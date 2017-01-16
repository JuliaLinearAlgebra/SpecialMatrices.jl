c = [1, 2 + 3.1im, 3, 4 - 2.5im]
C = Circulant(c)
x = convert(Vector{Complex{Float64}}, randn(size(C, 1)))
xt = transpose(x)
Cr = Circulant(randn(4))

# Important to check this as everything else depends upon it working correctly.
function check_full(C=C)
    Ĉ = [C.c circshift(C.c, 1) circshift(C.c, 2) circshift(C.c, 3)]
    all(C == Ĉ)
end
@test check_full()

# Check that, upon copying, the correct semantics are respected.
function check_copy(C=C)
    D = copy(C)
    D.c[1] = 5
    C.c[1] == 5
end
@test check_copy()

# Memory should not be shared between the result and the original.
function check_conj(C=C)
    D = conj(C)
    all(conj(full(C)) == full(conj(C))) && !all(full(D) == full(C))
end
@test check_conj()

# Memory should be shared between the new and the original.
function check_conj!(C=C)
    copy_C = deepcopy(C)
    D = conj!(copy_C)
    all(conj(full(C)) == full(D)) && all(full(D) == full(copy_C))
end
@test check_conj!()

function check_transpose(C=C)
    all(transpose(full(C)) == full(transpose(C)))
end
@test check_transpose()

function check_ctranspose()
    all(ctranspose(full(C)) == full(ctranspose(C)))
end
@test check_ctranspose()

function check_add(C=C)
    all((full(C) + full(C)) == full(C + C)) &&
        all((full(C) + 5) == full(C + 5)) &&
        all((5 + full(C)) == full(5 + C))
end
@test check_add()

function check_minus(C=C)
    all(full(-C) == -full(C)) &&
        all((full(C) - full(C)) == full(C - C)) &&
        all((full(C) - 5) == full(C - 5))
end
@test check_minus()

function check_eigen(C=C, x=x)
    eigen = eigfact(C)
    Γ = Diagonal(eigen.values)
    U = eigen.vectors
    abs(sum(Ac_mul_B(U, (Γ * (U * x))) - full(C) * x)) < 1e-9
end
@test check_eigen(Cr)

check_left_dot(C=C, x=x) = abs(sum(C * x - full(C) * x)) < 1e-9
@test check_left_dot()

check_right_dot(C=C, xt=xt) = abs(sum(xt * C - xt * full(C))) < 1e-9
@test check_right_dot()

check_backslash(C=C, x=x) = abs(sum(C * (C \ x) - x)) < 1e-9
@test check_backslash()

check_forwardslash(C=C, xt=xt) = abs(sum((xt / C) * C - xt)) < 1e-9
@test check_forwardslash()

function check_inplace_diag_mult()
    x = randn(100)
    d = randn(100)
    D = Diagonal(d)
    tmp = D * x
    all(tmp == A_mul_B!(D, x))
end
check_inplace_diag_mult()