import string
import sys
base_filename = sys.argv[1]

resolution = "0256"
initfile = open("%s_init_%s.txt" % (base_filename,resolution),"r").readlines()
fourthfile = open("%s_fourth_%s.txt" % (base_filename,resolution),"r").readlines()
midfile = open("%s_mid_%s.txt" % (base_filename,resolution),"r").readlines()
threefourthfile = open("%s_three_fourth_%s.txt" % (base_filename,resolution),"r").readlines()
endfile = open("%s_end_%s.txt" % (base_filename,resolution),"r").readlines()
solfile = open("%s_solution_%s.txt" % (base_filename,resolution),"r").readlines()


out = ["%s,%s,%s,%s,%s,%s,%s\n" % (string.split(i)[0],
                            string.split(i)[1],
							string.split(f)[1],
							string.split(m)[1],
							string.split(t)[1],
							string.split(e)[1],
							string.split(s)[1]) for (i,f,m,t,e,s) 
							in zip(initfile,fourthfile,midfile,threefourthfile,endfile,solfile) ]

outfile = open("%s_consolidated_%s.txt" % (base_filename,resolution),"w")
outfile.writelines(out)


