import socket, json

#
# Implement just enough asynchronous HTTP to let us issue Tuber calls.
# This does not use urllib or something because these libraries do not let
# you issue a request and get the reply later.
#

class TuberClient(object):
    def __init__(self, hostname, timeout=None):
        self.host = hostname
        self.socket = None
        self.timeout = timeout
    def OpenSocket(self):
        self.socket = socket.socket()
        self.socket.connect((self.host, 80))
        self.socket.settimeout(self.timeout)
    def ReadLineFromSocket(self):
        out = str()
        while not out.endswith('\n'):
            out += self.socket.recv(1).decode()
        return out
    def CallMethod(self, obj, method, wait_for_reply=True):
        if self.socket is None:
            self.OpenSocket()

        request = 'POST /tuber HTTP/1.1\r\n'
        request += 'Connection: keep-alive\r\n'
        request += 'Host: %s\r\n' % self.host
        request += 'Content-Type: application/json\r\n'

        if not isinstance(method, str) and hasattr(method, '__iter__'):
            jsonbits = '[%s]' % ', '.join(['{"object": "%s", "method": "%s"}' % (obj, m) for m in method])
        else:
            jsonbits = '[{"object": "%s", "method": "%s"}]' % (obj, method)
        request += 'Content-Length: %d\r\n\r\n' % len(jsonbits)
        request += jsonbits

        self.socket.send(request.encode())

        if wait_for_reply:
            return self.GetReply()
    def GetReply(self):
        toread = None
        chunked = False
        connection_atend = 'close'
        while True:
            header = self.ReadLineFromSocket().strip()
            if header == '':
                break
            header = header.split(':')
            if header[0] == 'Content-Length':
                toread = int(header[1].strip())
            if header[0] == 'Transfer-Encoding' and header[1].strip() == 'chunked':
                chunked = True
            if header[0] == 'Connection':
                connection_atend = header[1].strip()
        assert(toread is not None or chunked)
        assert(toread is None or chunked is False)
        out = ''

        if chunked:
            still_chunking = True
            while still_chunking:
                chunksize = self.ReadLineFromSocket().strip()
                toread = int(chunksize.split(';')[0], base=16)
                if toread == 0:
                    still_chunking = False
                while toread > 0:
                    chunk = self.socket.recv(toread).decode()
                    assert(len(chunk) > 0)
                    out += chunk
                    toread -= len(chunk)
                lf = self.ReadLineFromSocket().strip()
                assert(len(lf)) == 0
        else:
            while toread > 0:
                chunk = self.socket.recv(toread).decode()
                assert(len(chunk) > 0)
                out += chunk
                toread -= len(chunk)
        if connection_atend == 'close':
            self.socket.close()
            self.socket = None
        return json.loads(out)

