# Original Sources for MC (:= Medcouple) from the original author Guy Bris

As of May 2018, still available from
	https://wis.kuleuven.be/stat/robust/Programs/MEDCOUPLE/mc.zip

Where zipinfo shows :

Filemode    | Length |  Date      |     Time |      File
----------- | ----- :| ----------:| ------- :| ------------------
 -rw-rw-rw- | 21906 | 11-May-2002 | 09:15:50 | MC/spmc.c
 -rw-rw-rw- |  2949 | 25-Oct-2002 | 04:30:14 | MC/mcplus.ssc
 -rw-rw-rw- |   153 | 11-May-2002 | 10:14:28 | MC/mc.ssc
 -rw-rw-rw- |  5745 | 11-May-2002 | 09:39:38 | MC/spmc.obj
 -rw-rw-rw- |   795 | 16-Apr-2003 | 18:57:16 | MC/mc.txt
----------- | ----- | ----------- | -------- | -----------------
 -rw-rw-rw- |   892 |  3-Jul-2001 | 20:13:58 | MC/matlab/mc.c
 -rw-rw-rw- |  8192 |  2-Aug-2001 | 19:42:38 | MC/matlab/mc.dll
 -rw-rw-rw- | 21834 |  3-Jul-2001 | 20:14:04 | MC/matlab/mlmc.c
----------- | ----- | ----------- | -------- | -----------------
 -rw-rw-rw- |   106 |  5-Sep-2002 | 23:27:50 | MC/splus6/chapter.mif
 -rw-rw-rw- |  1226 |  5-Sep-2002 | 23:27:50 | MC/splus6/make.mak
 -rw-rw-rw- |   121 |  5-Sep-2002 | 23:28:04 | MC/splus6/S.def
 -rw-rw-rw- | 37888 |  5-Sep-2002 | 23:28:38 | MC/splus6/S.dll
----------- | ----- | ----------- | -------- | -----------------
            | 101807|             |          |   15 files
----------- | ----- | ----------- | -------- | -----------------

and I have removed the binary files `spmc.obj`, `mc.dll`, and `S.dll`
and kept the rest in this directory `robustbase/misc/MC/`

## Notes :

0. `mc.txt` gives the original overview by  Guy Bris (U Antwerpen, BE)
1. There's pure S (splus) version in `./mcplus.ssc`

2. The two main C (C++) source files are

    - `spmc.c`  {S-plus MC}
    - `mlmc.c`  {Matlab MC}

 and they very similar indeed, and contain very similar code
