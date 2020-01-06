using BinaryProvider
cd(@__DIR__) do
	@info("Testing BinaryProvider.unpack().")
	unpack("temp.gz", pwd())
	@assert(isfile("file 1.txt"))
	@assert(isfile("file 2.txt"))
end
