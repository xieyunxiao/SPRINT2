import sys,os

repeat_path=sys.argv[1]
low_com_path=sys.argv[2]
fi=open(repeat_path)
fo=open(low_com_path,'w')
while 1:
    lines = fi.readlines(10000)
    if not lines:
        break
    for line in lines:
        if not line:
            break

        this_line=line.strip("\n")
        if 'Simple_repeat' in this_line or 'Low_complexity' in this_line:
            fo.write(this_line+"\n")
fi.close()
fo.close()