H = Hilbert(10)
@test H[10,10] == 1//(10+10-1)
@test inv(H) == InverseHilbert(10)
