from __future__ import print_function
import copy
import os
import sys
if sys.version_info[0] >= 3:
        import urllib.parse as urlparse
        from urllib.error import HTTPError
        from urllib.request import urlopen, Request
else:
        import urlparse
        from urllib2 import HTTPError
        from urllib2 import urlopen, Request
import base64
import tempfile
import shutil
import subprocess
import re

from spt3g.core import G3Reader, G3Writer, G3FrameType
from spt3g.core import g3logging as logging

class GridFTPStager(object):
    """
    Provides to subprocess calls to stage files in and out using GridFTP
    
    Handles ftp:// and gsiftp:// URLs
    
    .. note:: GridFTP requires that you have a proxy certificate either in the
              standard location or in the location specified by the environment
              variable X509_USER_PROXY. See the `Globus Toolkit documentation 
              <http://toolkit.globus.org/toolkit/docs/4.1/admin/docbook/gtadmin-env-var.html#id2565277>`_
              for more information. You will also need to 
              `obtain a user certificate `_.
    """
    def __init__(self, globus_url_copy='globus-url-copy', options=['-nodcau', '-rst', '-cd']):
        super(type(self), self).__init__()
        self.globus_url_copy = globus_url_copy
        self.options = list(options)

    def transfer_file_in(self, url, local_path):
        """
        Transfers files to local directory from remote site
        
        Keyword arguments:
        url -- remote location of the file as a URL (no default)
        local_path -- path to local directory where file show be put (no default)
        """
        local_path = os.path.join(os.path.abspath(local_path), 
                                  os.path.basename(url))
        proc = subprocess.Popen([self.globus_url_copy] + self.options + [url, 
                                 'file://'+os.path.abspath(local_path)], 
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            logging.log_fatal("globus-url-copy failed with status %d: %s" % \
                              (proc.returncode, stderr.strip()), 
                              unit="GridFTPStager")
        else:
            logging.log_info("Download finished: %s to %s" % (url, local_path), 
                             unit="GridFTPStager")
        return local_path
    
    def transfer_file_out(self, local_path, url):
        """
        Transfers files to local directory from remote site
        
        Keyword arguments:
        local_path -- path to local directory where file is located (no default)
        url -- remote location of the file as a URL (no default)
        """
        proc = subprocess.Popen([self.globus_url_copy] + self.options + \
                                ['file://'+os.path.abspath(local_path), url], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            logging.log_fatal("globus-url-copy failed with status %d: %s" % \
                              (proc.returncode, stderr.strip()), 
                              unit="GridFTPStager")
        else:
            logging.log_info("Upload finished: %s to %s" % (local_path, url), 
                             unit="GridFTPStager")

class HTTPStager(object):
    """
    Provides to subprocess calls to stage files in and out using http. 
    
    Handles http://, https://, ftp://, and file:// URLs
    """
    def __init__(self, blocksize=2**16, ssl=None):
        self.blocksize = blocksize
        self.ssl = ssl if ssl else {}

    def strip_auth(self, url):
        """
        urlopen() doesn't support inline auth. Strip it out and
        construct the appropriate header by hand.
        
        Keyword arguments:
        url -- url with authenication tokens
        """
        parsed = urlparse.urlparse(url)
        auth = None
        if '@' in parsed.netloc:
            auth, netloc = parsed.netloc.split('@')
            auth = base64.encodestring(auth.encode('utf-8'))[:-1]
            parts = list(parsed)
            parts[1] = netloc
            url = urlparse.ParseResult(*parts).geturl()
        req = Request(url)
        if auth is not None:
            req.add_header('Authorization', 'Basic ' + auth.decode('utf-8'))
        return req

    def transfer_file_in(self, url, local_path):
        """
        Transfers files to local directory from remote site
        
        Keyword arguments:
        url -- remote location of the file as a URL (no default)
        local_path -- path to local directory where file show be put (no default)
        """
        # parse the URL
        # use "file" as the default scheme
        parsed_url = urlparse.urlparse(url, scheme="file")
        input_path = parsed_url.path
        output_path = os.path.join(local_path, os.path.basename(input_path))
        # copy the file
        if parsed_url.scheme == "file":
            if parsed_url.netloc != "":
                logging.log_fatal("%s is not a valid URI (file:// URIs must use absolute paths)" % url, 
                                  unit="HTTPStager")
            input_path = parsed_url.path
            logging.log_info("Copying file %s to %s" % (input_path, output_path),
                             unit="HTTPStager")
            shutil.copyfile(input_path, output_path)
            logging.log_info("File copied: %s to %s" % (input_path, output_path),
                             unit="HTTPStager")
        else:
            logging.log_info("Downloading %s to %s" % (url, output_path), 
                             unit="HTTPStager")
            f = None
            output_file = None
            try:
                output_file = open(output_path, "wb")
                try:
                    f = urlopen(self.strip_auth(url),**self.ssl)
                except TypeError:
                    f = urlopen(self.strip_auth(url))
                while True:
                    block = f.read(self.blocksize)
                    output_file.write(block)
                    if len(block) < self.blocksize:
                        break
                logging.log_info("Download finished: %s to %s" % (url, local_path), 
                                 unit="HTTPStager")
            except Exception as e:
                if os.path.exists(local_path):
                    os.remove(local_path)
                logging.log_fatal("Downloading %s: %s" % (url, str(e)), 
                                  unit="HTTPStager")
            finally:
                if f is not None:
                    f.close()
                if output_file is not None:
                    output_file.close()
        return output_path
    
    def transfer_file_out(self, local_path, url):
        """
        Transfers files to local directory from remote site. Only works for file:// handle. 
        Will we ever have an HTTP drop box? Else HTTP, etc. will not be supported
        
        Keyword arguments:
        local_path -- path to local directory where file is located (no default)
        url -- remote location of the file as a URL (no default)
        """
        parsed_url = urlparse.urlparse(url, scheme="file") # use "file" as the default scheme
        if parsed_url[0] not in ["http", "https", "ftp", "file"]:
            logging.log_fatal("Cannot handle URL scheme \"%s\": %s" % (parsed_url[0], url), 
                              unit="HTTPStager")
        # copy the file
        if parsed_url[0] == "file":
            output_path = parsed_url[2]
            if parsed_url.netloc != "":
                logging.log_fatal("%s is not a valid URI (file:// URIs must use absolute paths)" % url, 
                                  unit="HTTPStager")
            logging.log_info("Copying file %s to %s" % (local_path, output_path), 
                             unit="HTTPStager")
            shutil.copyfile(local_path, output_path)
            logging.log_info("File copied: %s to %s" % (local_path, output_path), 
                             unit="HTTPStager")
        else:
            logging.log_fatal("Can't upload to %s" % url.scheme)

class G3FileStager(object):
    """
    Staging mechanism class.
    
    Gets the input, output, and (optionally) temporary directory paths.
    Decides which transfer protocol to use depending on the file path.
    Provides function to stage files in and out. 
    """
    def __init__(self, input_path=None, output_path=None, 
                 tmp_dir=None, simulation=False, 
                 cleanup=True):
        if input_path is None:
            logging.log_fatal("Need to provide an input file")
        if output_path is None:
            logging.log_fatal("Need to provide an output file")
        if isinstance(input_path, list) or isinstance(input_path, str):
            self.input_path = input_path
        else:
            raise RuntimeError("Please provided the input file(s) as a list or string")
        if isinstance(output_path, list):
            self.output_path = output_path
        elif isinstance(output_path, str):
            self.output_path = [output_path]
        else:
            raise RuntimeError("Please provided the output file(s) as a list or string")
        self.tmp_dir = os.path.join(self.get_local_scratch_dir(tmp_dir), "input")
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        self.cleanup=cleanup
        self.local_input_files = []

    
    def stage_files_in(self):
        if isinstance(self.input_path, list):
            local_files = [self.stage_file_in(f) for f in self.input_path]
        else:
            local_files = [self.stage_file_in(self.input_path)]
        return local_files
        
    def stage_file_in(self, input_file):
        return self.determine_protocol(input_file).transfer_file_in(input_file, self.tmp_dir)
        
    def stage_file_out(self, local_output_filename, remote_output_filename):
        self.determine_protocol(remote_output_filename).transfer_file_out(local_output_filename, remote_output_filename)
        if self.cleanup:
            os.remove(local_output_filename)
        
    def get_local_g3_filenames(self):
        if len(self.local_input_files) == 0: 
            self.local_input_files = self.stage_files_in()
        return [f for f in self.local_input_files if ".g3" in os.path.basename(f)]
        
    def get_remote_g3_output_filename(self):
        remote_outfiles = [f for f in self.output_path if ".g3" in os.path.basename(f)]
        if len(remote_outfiles) > 1:
            logging.log_warn("Ignoring multiple output .g3 files. Using %s as the output file name" % remote_outfiles[0],
                             unit="G3FileStager")
        if len(remote_outfiles) < 1:
            logging.log_fatal("Please provide a .g3 output file.", 
                              unit="G3FileStager")
        return remote_outfiles[0]
    
    def get_local_output_filename(self, filename):
        return os.path.join(self.tmp_dir, os.path.basename(filename))
        
    def get_local_nong3_filenames(self):
        if len(self.local_input_files) == 0: 
            self.local_input_files = self.stage_files_in()
        return [f for f in self.local_input_files if ".g3" not in os.path.basename(f)]
        
    def determine_protocol(self, filename):
        protocol_id = urlparse.urlparse(filename).scheme
        if protocol_id in ["gsiftp", "ftp"]:
            return GridFTPStager()
        elif protocol_id in ["http", "https", "ftp", "file"]:
            return HTTPStager()
        else:
            raise NotImplementedError("The %s protocal is not supported" % protocol_id)
    
    def try_to_make_scratch_dir(self, basename, fullname):
        if not os.path.isdir(basename):
            return False
        if os.path.isdir(fullname):
            return True
        try:
            os.mkdir(fullname)
        except OSError:
            return False
        return True

    def get_local_scratch_dir(self, tmp_dir):
        if tmp_dir is not None:
            return tmp_dir
        # find a local staging directory
        elif "_CONDOR_SCRATCH_DIR" in os.environ:
            # works on condor (especially npx4)
            staging_directory = os.environ["_CONDOR_SCRATCH_DIR"]
        elif "TMPDIR" in os.environ:
            # works on some PBS/TORQUE nodes
            staging_directory = os.environ["TMPDIR"]
        else:
            # try some known locations on interactive nodes
            import pwd
            current_username = pwd.getpwuid(os.getuid()).pw_name
            if self.try_to_make_scratch_dir('/scratch', '/scratch/'+current_username):
                staging_directory = '/scratch/'+current_username
            elif self.try_to_make_scratch_dir('/global/scratch', '/global/scratch/'+current_username):
                staging_directory = '/global/scratch/'+current_username
            else:
                staging_directory = os.getcwd()
                logging.log_info("Cannot find a suitable scratch directory on this machine; falling back to the current working directory (%s). If this is not what you want, set a different path with dataio.set_local_scratch_dir(path)." % staging_directory, unit="G3Reader")
            return staging_directory
    
@core.indexmod
class G3URLReader(object):
    """
    G3Reader that takes a G3FileStager() instead of filenames
    """
    def __init__(self, stager = None):
        self.infiles = stager.get_local_g3_filenames()
        print(self.infiles)
        self.reader = G3Reader(self.infiles.pop(0))
    
    def __call__(self, frame):
        assert(frame is None) # Needed for driving module
        outframe = self.reader(None)
        if len(outframe) == 0:
            while len(self.infiles) > 0 and len(outframe) == 0: 
                self.reader = G3Reader(self.infiles.pop(0))
                outframe = self.reader(None)
            return outframe
        else:
            return outframe

@core.indexmod
class G3URLWriter(object):
    """
    G3Writer that takes a G3FileStager() instead of filenames
    """
    def __init__(self, stager = None):
        self.stager = stager
        self.remote_outfile = stager.get_remote_g3_output_filename()
        self.local_outfile = stager.get_local_output_filename(self.remote_outfile)
        self.writer = G3Writer(self.local_outfile)
    
    def __call__(self, frame):
        if frame.type != G3FrameType.EndProcessing:
            self.writer(frame)
        else:
            del self.writer
            self.stager.stage_files_out(self.local_outfile, self.remote_outfile)
            
