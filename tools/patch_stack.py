import struct,sys
src,dst,mb=sys.argv[1],sys.argv[2],int(sys.argv[3])
b=bytearray(open(src,'rb').read())
pe=struct.unpack_from('<I',b,0x3c)[0]
opt=pe+24
magic=struct.unpack_from('<H',b,opt)[0]
assert magic==0x20b
off=opt+72  # SizeOfStackReserve (PE32+)
print('old reserve',struct.unpack_from('<Q',b,off)[0])
struct.pack_into('<Q',b,off,mb*1024*1024)
open(dst,'wb').write(b)
