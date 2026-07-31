function aplicaCDC_u(Ht, Gt, dad)
    n = size(dad.NOS, 1)
    tipo = tipoCDC(dad)
    valoresconhecidos = valorCDC(dad)

    normal_fonte = dad.normal
    CDCuc = tipo[1:n] .== 0
    F, Fx, Fy = BEM.montaFs(dad.NOS[CDCuc, :], [dad.NOS; dad.pontos_internos])

    Gu =
        -dad.k *
        Gt[:, CDCuc] *
        (Fx .* normal_fonte[CDCuc, 1] + Fy .* normal_fonte[CDCuc, 2])
    # @infiltrate
    A = Ht - Gu
    b =
        Gt[:, (1:n)[tipo[1:n]]] * valoresconhecidos[tipo[1:n]] -
        Ht[:, (1:n)[.!tipo[1:n]]] * valoresconhecidos[.!tipo[1:n]]

    tipouc = tipo .== 0

    Add = A[tipo, tipo]
    Adc = A[tipo, tipouc]
    bd = b[tipo]
    ud=Add \ (bd-Adc*valoresconhecidos[tipo[1:n] .== 0])
    T=similar(b)
    q=similar(valoresconhecidos)
    T[tipouc] = valoresconhecidos[tipo[1:n] .== 0]
    T[tipo] = ud

    q[tipouc[1:n]] =
        -dad.k * (Fx .* normal_fonte[CDCuc, 1] + Fy .* normal_fonte[CDCuc, 2]) * T
    q[tipo[1:n]] = valoresconhecidos[tipo[1:n]]
    T, q
end

function aplicaCDC_alpha(Ht, Gt, dad)
    n = size(dad.NOS, 1)
    tipo = tipoCDC(dad)
    valoresconhecidos = valorCDC(dad)

    normal_fonte = dad.normal
    CDCuc = tipo[1:n] .== 0
    F, dFdx, dFdy, ~ =
        BEM.montaF([dad.NOS; dad.pontos_internos], [dad.NOS; dad.pontos_internos])

    dFdn = -(dFdx[:, 1:n]' .* normal_fonte[:, 1] + dFdy[:, 1:n]' .* normal_fonte[:, 2])

    Ga = -dad.k * Gt[:, CDCuc] * dFdn[CDCuc, :]
    Ha = Ht[:, tipo] * F[tipo, :]
    A = Ha - Ga
    b = Gt[:, .!CDCuc] * valoresconhecidos[.!CDCuc]-Ht[:, .!tipo] * valoresconhecidos[CDCuc]


    a = A \ b
    T = F * a
    q = -dad.k * dFdn * a
    T, q
end


function corrige_autovalor(G, H, M, dad)
    n = nc(dad) + ni(dad)
    ind_ud = trues(n)
    ind_qd = falses(nc(dad))
    for elem_i in dad.ELEM  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[1] == 0
            ind_qd[ind_elem] .= 1
            ind_ud[ind_elem] .= 0
        end
    end
    guu = G[(1:nc(dad))[ind_qd], ind_qd]
    gqu = G[ind_ud, ind_qd]
    huq = H[(1:nc(dad))[ind_qd], ind_ud]
    hqq = H[ind_ud, ind_ud]
    muq = M[(1:nc(dad))[ind_qd], ind_ud]
    mqq = M[ind_ud, ind_ud]
    Hb = hqq - gqu * (guu \ huq)
    Mb = mqq - gqu * (guu \ muq)
    a, v = eigen(Hb, -Mb)
    v_wn = sort(sqrt.(abs.(real(a))))
end


function corrige_autovalor_u(Gt, Ht, M, dad)
    n = size(dad.NOS, 1)

    ud = trues(n + ni(dad))
    qd = falses(nc(dad))

    for elem_i in dad.ELEM  # Laço dos pontos fontes
        ind_elem = elem_i.indices
        if elem_i.tipoCDC[1] == 0
            ud[ind_elem] .= 0
            qd[ind_elem] .= 1
        end
    end

    normal_fonte = dad.normal
    F, Fx, Fy = BEM.montaFs(dad.NOS[qd, :], [dad.NOS; dad.pontos_internos])

    Gu = -dad.k * Gt[:, qd] * (Fx .* normal_fonte[qd, 1] + Fy .* normal_fonte[qd, 2])
    # @infiltrate
    Hu = Ht - Gu
    a, v = eigen(Hu[ud, ud], -M[ud, ud])
    v_wn = sort(sqrt.(abs.(real(a))))
end
