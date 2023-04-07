from multiprocessing import get_context
import socket, pickle, errno, struct, time

from spt3g.core import G3FrameType, G3Frame

# Require fork to avoid pickling errors
ctx = get_context("fork")

class Subproc(ctx.Process):
    '''
    Run a module in a subprocess, using python multiprocessing to proxy
    frames to it. If more than maxqueuelen frames are queued on the
    remote module, will block until there is space.
    '''
    def __init__(self, target=None, name=None, maxqueuelen=50):
        '''
        Set up a multiprocessing shim for the pipeline module <target>
        under the process name <name> (can be None).
        '''
        ctx.Process.__init__(self, name=name)
        self.targetmod = target

        self.queue = socket.socketpair()
        self.queue[0].setblocking(False)
        self.callsqueued = 0 # Count of frames outstanding at child

        self.maxqueued = maxqueuelen
        self.started = False
        self.outbuf = b''
    def __call__(self, frame):
        # Processing flow:
        # - Frames added to input queue (this will alert the remote
        #   process).
        # - Check output queue from the remote process and inject
        #   any frames there into the pipeline.

        if not self.started:
            self.started = True
            self.start()

        outbuf = frame.__getstate__()[1]
        self.outbuf += struct.pack('i', len(outbuf)) + outbuf
        self.callsqueued += 1

        # Loop on output queue so long as
        # (a) we are expecting data back and
        # (b) (i) we need to finalize processing (an EndProcessing frame
        #     is here), (ii) we need to drain the queue, *or* (iii)
        #     there are pending frames from the remote end
        #
        # In general, attempt to make as much forward progress on both
        # sending and receiving data as possible (or necessary) in the
        # current moment. If we can't make progress, and don't need to,
        # save state and try again next time.

        allout = []
        inbuf = b''
        icount = 0
        while icount < 10 or \
         (frame.type == G3FrameType.EndProcessing and self.callsqueued > 0) or \
         self.callsqueued >= self.maxqueued:
            if len(self.outbuf) > 0:
                try:
                    resid = self.queue[0].send(self.outbuf)
                except socket.error as e:
                    # Just continue on EAGAIN
                    if e.errno != errno.EAGAIN:
                        raise
                else:
                    self.outbuf = self.outbuf[resid:]

            if self.callsqueued > 0:
                try:
                    inbuf = self.queue[0].recv(4)
                except socket.error as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    # Fall through and do nothing if queue empty
                else:
                    # Unqueue a packet (length + data) from remote end
                    bytestoread = struct.unpack('i', inbuf)[0]
                    inbuf = b''
                    self.queue[0].setblocking(True)
                    while bytestoread > len(inbuf):
                        inbuf += self.queue[0].recv(bytestoread - len(inbuf))
                    self.queue[0].setblocking(False)
                                
                    allout += pickle.loads(inbuf)
                    self.callsqueued -= 1
                    icount = 0 # Forward progress! Try again quick

            icount += 1

        if frame.type == G3FrameType.EndProcessing:
            self.join()
        return allout
    def run(self):
        # Processing flow in worker process:
        # - Receive data from the remote end (blocking)
        # - Process data
        # - Return results, packaged as a list of frames
        #
        # Note that this must return *something* in response to every
        # frame to keep the outstanding frame accounting in the master
        # process working (an empty list is fine).

        while True:
            framelen = struct.unpack('i', self.queue[1].recv(4))[0]
            framebuf = b''
            while len(framebuf) < framelen:
                framebuf += self.queue[1].recv(framelen - len(framebuf))
            frame = G3Frame()
            frame.__setstate__(({}, framebuf))

            outframes = self.targetmod(frame)

            if outframes is None:
                out = [frame]
            elif isinstance(outframes, G3Frame):
                out = [outframes]
            elif hasattr(outframes, '__iter__'):
                out = list(outframes)
            elif outframes:
                out = [frame]
            else:
                out = []

            outbuf = pickle.dumps(out)
            outbuf = struct.pack('i', len(outbuf)) + outbuf
            self.queue[1].send(outbuf)
                
            if frame.type == G3FrameType.EndProcessing:
                return

