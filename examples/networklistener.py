from spt3g import core, dfmux, networkstreamer, hk

pipe = core.G3Pipeline()
pipe.Add(networkstreamer.G3NetworkReceiver, 
         hostname = 'laphroaig.berkeley.edu', 
         port = 8675)



pipe.Add(core.Dump)

pipe.Run()
