from spt3g.core import G3Reader

class G3File(object): 
	'''Iterable class for G3 files, as created by G3Writer. Loop through frames by doing something like:
	f = core.G3File('/path/to/file.g3')
	for frame in f:
		print( frame )

	An entire file can also be read into an indexable list by doing:
	f = list(core.G3File('/path/to/file.g3'))
	'''
	def __init__(self, path):
		self.reader = G3Reader(path)
	def __iter__(self):
		return self
	def next(self):
		frames = self.reader.Process(None)
		if len(frames) == 0:
			raise StopIteration('No more frames in file')
		if len(frames) > 1:
			raise ValueError('Too many frames returned by reader')
		return frames[0]
	__next__ = next

