c = [1, 2.5 + 3.1im, 3, 4 - 2.5im]
C = Circulant(c)
x = convert(Vector{Complex{Float64}}, randn(size(C, 1)))
xt = transpose(x)
Cr = Circulant(randn(4))

ϵ = 1e-12

# Important to check this as everything else depends upon it working correctly.
function check_full(C=C)
    Ĉ = [C.c circshift(C.c, 1) circshift(C.c, 2) circshift(C.c, 3)]
    all(full(C) == Ĉ)
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
    abs(sum(Ac_mul_B(U, (Γ * (U * x))) - full(C) * x)) < ϵ
end
@test check_eigen(Cr)

function check_to_circ(C=C)
    out = tocirc(eigfact(C))
    abs(sum(full(C) - full(out))) < ϵ
end
@test check_to_circ()

check_left_dot(C=C, x=x) = abs(sum(C * x - full(C) * x)) < ϵ
@test check_left_dot()

check_right_dot(C=C, xt=xt) = abs(sum(xt * C - xt * full(C))) < ϵ
@test check_right_dot()

check_backslash(C=C, x=x) = abs(sum(C * (C \ x) - x)) < ϵ
@test check_backslash()

check_forwardslash(C=C, xt=xt) = abs(sum((xt / C) * C - xt)) < ϵ
@test check_forwardslash()

function check_left_dot!(C=C, x=x)
    x_loc = deepcopy(x)
    x_out = deepcopy(x_loc)
    tmp = C * x_loc
    A_mul_B!(x_out, C, x_loc)
    abs(sum(x_out - tmp)) < ϵ
end
@test check_left_dot!()

function check_right_dot!(C=C, xt=xt)
    xt_loc = deepcopy(xt)
    tmp = xt_loc * C
    out = A_mul_B!(xt_loc, xt_loc, C)
    abs(sum(x - tmp)) < ϵ
end

check_inv(C=C) = abs(sum(inv(full(C)) - full(inv(C)))) < ϵ
@test check_inv()

check_det(C=C) = abs(det(full(C)) - det(C)) < ϵ
@test check_det()

check_logdet(C=C) = abs(logdet(full(C)) - logdet(C)) < ϵ
@test check_det()

check_real(C=C) = abs(sum(real(full(C)) - full(real(C)))) < ϵ
@test check_real()

function check_round(C=C)
    abs(sum(round(full(real(C))) - full(round(real(C))))) < ϵ
end
@test check_round()

function check_convert(C=C)
    from_full = convert(Matrix{Int64}, full(round(real(C))))
    from_circ = full(convert(Circulant{Int64}, real(round(C))))
    abs(sum(from_full - from_circ)) < ϵ
end
@test check_convert()
