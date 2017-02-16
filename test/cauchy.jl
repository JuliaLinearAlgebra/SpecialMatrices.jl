@test maxabs(1 ./ Cauchy(3) - [2 3 4; 3 4 5; 4 5 6]) < 1e-14
@test maxabs(1 ./ full(Cauchy(3)) - [2 3 4; 3 4 5; 4 5 6]) < 1e-14
@test maxabs(1 ./ Cauchy(1:3, [2im, 10]) - [1+2im 11; 2+2im 12; 3+2im 13]) < 1e-14
