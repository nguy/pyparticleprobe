__declspec(dllexport) long decompress_2d(PUCHAR bytes, long bytecount, PUCHAR uncompressed_bytes, long max_uncompressed_bytes, long *position)
{
long cur_in = 0;
long cur_out = 0;
unsigned char a_byte;
unsigned char iterator = 0;
long idx;
   long bad_msg = 0;
long use_a_byte = 0;
   long uncompressed_bytecount = 0;
   // The first byte is always the index of the first RLE header byte
cur_in = 0;
while (cur_in < bytecount) {
    if (cur_in > bytecount) {
bad_msg = 1;
        break;  // Want out of While
    }
    use_a_byte = 1;
           if ((bytes[cur_in] & 0xFF) & 0x20) {
                  cur_in++;
continue;
}
           else if ((bytes[cur_in] & 0xFF) > 0x7F) {
        a_byte = 0x00;
        //iterator = bytes[cur_in] - 128;
        iterator = bytes[cur_in] & 0x7F;
cur_in++;
}
else if ((bytes[cur_in] & 0xFF) > 0x3F) { // > 63
        a_byte = 0xFF;
//iterator = (bytes[cur_in] - 64) & 0xFF;
        iterator = bytes[cur_in] & 0xBF;
                  cur_in++;
}
    else {
        iterator = bytes[cur_in] & 0xFF;
        use_a_byte = 0;
        cur_in++;
        }
// > 127
if ((iterator & 0xFF) > 0x1F) { // > 31
//NiceMsgBox "Decompress:Iterator > 31 at " & cur_in - 1
Form DOC-0201 Rev B
© 2011 DROPLET MEASUREMENT TECHNOLOGIES, INC.
2 5
￼
Image Probe Data Reference Manual
￼}
              *position = cur_in;
    return (-31);  // Exit Sub
              }
if ((cur_out + iterator) >= max_uncompressed_bytes)
{ }
*position = cur_in;
return (uncompressed_bytecount);
        //For idx = iterator To 0 Step -1
        for (idx = iterator ; idx >= 0 ; idx--)
                      {
                      if (cur_out >= max_uncompressed_bytes) {
                             // uncompressed_bytes buffer is full
                             uncompressed_bytecount = cur_out - 1;
                             *position = cur_in;
                             return (uncompressed_bytecount);
                      }
            if (use_a_byte) {
                uncompressed_bytes[cur_out] = a_byte;
                cur_out++;
} else {
                if (cur_in > bytecount) {
                    bad_msg = 1;
break; }
                uncompressed_bytes[cur_out] = bytes[cur_in];
                cur_in++;
                cur_out++;
}
}
    if (bad_msg) {
               *position = cur_in;
        return (-1);
       }
    uncompressed_bytecount = cur_out - 1;
    *position = cur_in;
    return (uncompressed_bytecount);
}