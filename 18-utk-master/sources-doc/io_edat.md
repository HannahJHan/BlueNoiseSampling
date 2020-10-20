# File format .edat


## Files


```
src/io/fileIO.hpp  
src/io/fileIO_model.hpp  
src/io/fileIO_edat.hpp
```


## Description


UnicornTK relies on the extension of the input file to determine its type. When the extension is `.edat`, the pointsets are stored in ASCII mode, with a 3 line header containing the extend of the pointset domain and the toricity of the pointset. The first line is the minimal point of the domain, the second line is the maximal point and the third line contains a boolean indicating is the pointset is toric or not It separates pointsets using "#".

## Examples


Pointset 2D, Toric, in domain $[0,1]^3$, with integer coordinates:  

```
0       0
1024    2187
1
0       0
512     729
256     1458
768     243
128     972
640     1701
384     486
896     1215
64      1944
576     81
320     810
832     1539
192     324
704     1053
448     1782
960     567
32      1296
544     2025
288     162
800     891
160     1620
672     405
416     1134
928     1863
96      648
```


Several pointsets 3D, Toric, in domain $[0,1]^3$, with floating point coordinates:

```
0.000000000000000       0.000000000000000       0.000000000000000
1.000000000000000       1.000000000000000       1.000000000000000
1
0.050395372335425       0.605445716671051       0.404234371009337
0.825739547336879       0.019014228372605       0.571619643962853
0.251084789226550       0.796739092282192       0.662392276959119
0.076733451587319       0.966740387484138       0.104033727921553
0.022851851774112       0.541269098301714       0.274100657150638
0.836165956649037       0.710473378354328       0.218262058035066
0.734431378481544       0.045624867092265       0.108525469691615
0.424232415953705       0.622989150374892       0.445836456596348
0.459817434535319       0.149720273248321       0.263248428763647
0.492075960732083       0.549053552688345       0.021870193542417
0.331534335234094       0.343634019220163       0.221821374311536
0.398141971416339       0.458952383943590       0.978517358419140
0.269279661418237       0.511562038332122       0.016162883494087
0.258357287444723       0.187212941280383       0.399687562696563
0.796217287843166       0.147599688299838       0.400614753225961
0.748264057735974       0.041631821552671       0.061114694238900
#
0.000000000000000       0.000000000000000       0.000000000000000
1.000000000000000       1.000000000000000       1.000000000000000
1
0.728055556474266       0.560123025104432       0.385942792609786
0.199736131401671       0.196849285437923       0.846230597897952
0.965195278163346       0.820035925093115       0.776832124632046
0.207031543414814       0.707106014878281       0.422963624918592
0.437438748645931       0.934516076914621       0.731540281030242
0.884149564868805       0.817359755006004       0.868729761072604
0.340246866303554       0.478473960766213       0.475246819778170
0.919041518335939       0.328492126969735       0.889999115348328
0.629857444816701       0.050157970760520       0.017919836337194
0.936117153134224       0.322178287972272       0.747246075828384
0.260712328427637       0.595219216448073       0.609722221644996
0.703291282920002       0.980232460419701       0.416594642544304
0.947554823464610       0.746043482969060       0.298958397540021
0.861714192401086       0.825613557320464       0.466442995114821
0.609747354781662       0.819148976779345       0.878929865983997
0.334699008226092       0.987974778280588       0.963871809412695
```