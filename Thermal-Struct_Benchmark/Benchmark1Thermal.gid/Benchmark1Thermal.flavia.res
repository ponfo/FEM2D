GiD Post Result File 1.0
 Result "Temperature" "Benchmark1Thermal"            1  Scalar onNodes
 Values
           1   25.0000000000000     
           2   25.0000000000000     
           3   25.0000000000000     
           4   25.0000000000000     
           5   25.0000000000000     
           6   25.0000000000000     
           7   25.0000000000000     
           8   25.0000000000000     
           9   25.0000000000000     
          10   25.0000000000000     
          11   25.0000000000000     
          12   25.0000000000000     
          13   25.0000000000000     
          14   25.0000000000000     
          15   25.0000000000000     
          16   25.0000000000000     
          17   25.0000000000000     
          18   25.0000000000000     
          19   25.0000000000000     
          20   25.0000000000000     
          21   25.0000000000000     
          22   25.0000000000000     
          23   25.0000000000000     
          24   25.0000000000000     
          25   25.0000000000000     
          26   25.0000000000000     
          27   25.0000000000000     
          28   25.0000000000000     
          29   25.0000000000000     
          30   25.0000000000000     
          31   25.0000000000000     
          32   25.0000000000000     
          33   25.0000000000000     
          34   25.0000000000000     
          35   25.0000000000000     
          36   25.0000000000000     
          37   25.0000000000000     
          38   25.0000000000000     
          39   25.0000000000000     
 End Values

GaussPoints "PointsFluxOnTriangs" ElemType Triangle
Number of GaussPoints: 4
Natural Coordinates: Given
        0.3333333333333333          0.3333333333333333
        0.2000000000000000          0.2000000000000000
        0.6000000000000000          0.2000000000000000
        0.2000000000000000          0.6000000000000000
End gausspoints
Result "FluxOnTriangs" "Benchmark1Thermal" 1 Vector onGaussPoints "PointsFluxOnTriangs"
Values
1          0.0000000000000039          0.0000000000000000
              0.0000000000000036         -0.0000000000000107
              0.0000000000000071         -0.0000000000000107
              0.0000000000000107         -0.0000000000000142
2         -0.0000000000000025          0.0000000000000000
             -0.0000000000000036          0.0000000000000000
              0.0000000000000030          0.0000000000000142
             -0.0000000000000142          0.0000000000000071
3         -0.0000000000000036          0.0000000000000000
             -0.0000000000000018         -0.0000000000000071
             -0.0000000000000053         -0.0000000000000071
             -0.0000000000000142         -0.0000000000000071
4         -0.0000000000000023          0.0000000000000071
              0.0000000000000000          0.0000000000000071
              0.0000000000000027          0.0000000000000036
             -0.0000000000000107          0.0000000000000142
5          0.0000000000000000          0.0000000000000000
              0.0000000000000000          0.0000000000000071
              0.0000000000000089          0.0000000000000036
              0.0000000000000000          0.0000000000000036
6          0.0000000000000032          0.0000000000000000
              0.0000000000000107         -0.0000000000000107
              0.0000000000000142         -0.0000000000000107
              0.0000000000000071         -0.0000000000000213
7          0.0000000000000000          0.0000000000000107
              0.0000000000000036          0.0000000000000213
              0.0000000000000089          0.0000000000000142
              0.0000000000000000          0.0000000000000000
8         -0.0000000000000010          0.0000000000000071
              0.0000000000000000         -0.0000000000000036
             -0.0000000000000071          0.0000000000000107
              0.0000000000000000          0.0000000000000071
9         -0.0000000000000036         -0.0000000000000107
             -0.0000000000000018         -0.0000000000000142
             -0.0000000000000053         -0.0000000000000071
              0.0000000000000000          0.0000000000000000
10          0.0000000000000051          0.0000000000000000
              0.0000000000000107          0.0000000000000071
              0.0000000000000096         -0.0000000000000107
             -0.0000000000000036          0.0000000000000000
11          0.0000000000000000          0.0000000000000071
              0.0000000000000036          0.0000000000000000
              0.0000000000000071         -0.0000000000000178
              0.0000000000000000          0.0000000000000036
12          0.0000000000000060          0.0000000000000000
              0.0000000000000036         -0.0000000000000107
              0.0000000000000071         -0.0000000000000107
              0.0000000000000071         -0.0000000000000071
End Values