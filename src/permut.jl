###########################
### Ranndom permutation ###
###########################

function permuteX(wlt::Array{Float64,2},taut::Array{Float64,4},t::Array{Float64,4},rhot::Float64,
                        x::Array{Float64,2},p::Array{Float64,2},taup::Array{Float64,2},
                        fty::Array{Float64,2},rhop::Float64,taus::Array{Float64,2},s::Array{Float64,2},rhos::Float64,
                        spatialwlt,nspat::Int,psfcbe::Array{Complex64,2})

    if spatialwlt[end] == "dirac"
        wlt = myidwt_dirac(wlt, nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
        b = fty + taup + rhop*p + taus + rhos*s
        wlt_b = wlt + b
        x = real(ifft(psfcbe.*fft(wlt_b)))
    else
        wlt = myidwt(wlt, nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
        b = fty + taup + rhop*p + taus + rhos*s
        wlt_b = wlt + b
        x = real(ifft(psfcbe.*fft(wlt_b)))
    end
    return x
end







function permuteP(x::Array{Float64,2},p::Array{Float64,2},taup::Array{Float64,2},
                  rhop::Float64,mask::Array{Float64,2})

        tmp = x-taup/rhop
        p = max(0,tmp).*mask

    return p

end


function permuteT(taut::Array{Float64,4}, t::Array{Float64,4}, rhot::Float64, x::Array{Float64,2},
                  spatialwlt, μt::Float64, nspat::Int)

    if spatialwlt[end] == "dirac"
        for b in 1:nspat-1
                    hx = dwt(x,wavelet(spatialwlt[b]))
                    tmp = hx - taut[:,:,1,b]/rhot
                    t[:,:,1,b] = prox_u(tmp,μt/rhot)
        end
                    hx = x
                    tmp = hx - taut[:,:,1,nspat]/rhot
                    t[:,:,1,nspat] = prox_u(tmp,μt/rhot)

    else
        for b in 1:nspat
                    hx = dwt(x,wavelet(spatialwlt[b]))
                    tmp = hx - taut[:,:,1,b]/rhot
                    t[:,:,1,b] = prox_u(tmp,μt/rhot)
        end
    end

    return t

end
