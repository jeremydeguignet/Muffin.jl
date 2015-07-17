# function parallelmuffin(wlt,taut,t,rhot,x,psf,psfadj,p,taup,fty,rhop,taus,s,rhos,mu,spatialwlt,μt,nspat)
#
#     wlt = myidwt(wlt, nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
#     b = fty + taup + rhop*p + taus + rhos*s
#     wlt_b = wlt + b
#
#     nxy = (size(x))[1]
#     nxypsf = (size(psf))[1]
#     psfcbe = zeros(Complex64,nxy,nxy)
#     psfpad = zeros(Float64,nxy,nxy)
#     psfpad[1:nxypsf,1:nxypsf] = psf[:,:]
#     psfcbe = 1./(abs(fft(psfpad)).^2+mu)
#     x = real(ifft(psfcbe.*fft(wlt_b)))
#
#     # tmp1 = 0.0
#     # tmp2 = zeros(Float64,nxy,nxy)
#     for b in 1:nspat
#                 hx = dwt(x,wavelet(spatialwlt[b]))
#                 tmp = hx - taut[:,:,1,b]/rhot
#                 t[:,:,1,b] = prox_u(tmp,μt/rhot)
#                 taut[:,:,1,b] = taut[:,:,1,b] + rhot*(t[:,:,1,b]-hx)
#                 # tmp1 = vecnorm([tmp2 (hx-(admmst.t)[:,:,z,b])],2)
#                 # tmp2 = (hx-(admmst.t)[:,:,z,b])
#     end
#     # tmp2[:] = 0
#
#
#     tmp = x-taup/rhop
#     p = max(0,tmp)
#     taup = taup + rhop*(p-x)
#
#     return wlt,x,t,taut,p,taup
#
# end

function parallelmuffin(wlt,taut,t,rhot,x,psf,psfadj,p,taup,fty,rhop,taus,s,rhos,mu,spatialwlt,μt,nspat)

    wlt = convert(Array,wlt)
    x = convert(Array,x)
    p = convert(Array,p)
    taup = convert(Array,taup)
    t = convert(Array,taut)
    taut = convert(Array,taut)

    wlt = myidwt(wlt, nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
    b = fty + taup + rhop*p + taus + rhos*s
    wlt_b = wlt + b

    nxy = (size(x))[1]
    nxypsf = (size(psf))[1]
    psfcbe = zeros(Complex64,nxy,nxy)
    psfpad = zeros(Float64,nxy,nxy)
    psfpad[1:nxypsf,1:nxypsf] = psf[:,:]
    psfcbe = 1./(abs(fft(psfpad)).^2+mu)
    x = real(ifft(psfcbe.*fft(wlt_b)))

    # tmp1 = 0.0
    # tmp2 = zeros(Float64,nxy,nxy)
    for b in 1:nspat
                hx = dwt(x,wavelet(spatialwlt[b]))
                tmp = hx - taut[:,:,1,b]/rhot
                t[:,:,1,b] = prox_u(tmp,μt/rhot)
                taut[:,:,1,b] = taut[:,:,1,b] + rhot*(t[:,:,1,b]-hx)
                # tmp1 = vecnorm([tmp2 (hx-(admmst.t)[:,:,z,b])],2)
                # tmp2 = (hx-(admmst.t)[:,:,z,b])
    end
    # tmp2[:] = 0


    tmp = x-taup/rhop
    p = max(0,tmp)
    taup = taup + rhop*(p-x)

    return wlt,x,t,taut,p,taup

end


@parallel for i in 2:nfreq+1
    z = i-1
    admmst.wlt[:,:,z],admmst.x[:,:,z],admmst.t[:,:,z,:],admmst.taut[:,:,z,:],admmst.p[:,:,z],admmst.taup[:,:,z] =
                                @fetchfrom(i,parallelmuffin(admmst.wlt[:,:,z], admmst.taut[:,:,z,:], admmst.t[:,:,z,:], rhot, admmst.x[:,:,z],
                                psfst.mypsf[:,:,z], psfst.mypsfadj[:,:,z], admmst.p[:,:,z], admmst.taup[:,:,z],
                                fty[:,:,z], rhop, admmst.taus[:,:,z], admmst.s[:,:,z], rhos, admmst.mu, spatialwlt,
                                μt, nspat))
end




image = randn(2048,2048,16);
psf = randn(2048,2048,16);
resultats = SharedArray(Float64,2048,2048,16);


@time @parallel for z in 1:16
    resultats[:,:,z] = (image[:,:,z].*psf[:,:,z])
    psf[:,:,z] = psf[:,:,z]*2
    println(sum(resultats[:,:,z]))
    # rr = RemoteRef(myid())
    # put!(rr,psf[:,:,z])
end

@time @parallel for z in 1:16
    resultats = @fetchfrom(z+1, gc())
    image = @fetchfrom(z+1, gc())
    psf = @fetchfrom(z+1, gc())
end














MAP = { (admmst.taut[:,:,z,:], admmst.t[:,:,z,:], rhot, admmst.x[:,:,z], psfst.mypsf[:,:,z],
         admmst.p[:,:,z], admmst.taup[:,:,z], admmst.fty[:,:,z], rhop, admmst.taus[:,:,z],
         admmst.s[:,:,z], rhos, mu, spatialwlt, μt, nspat) for z in 1:nfreq};


@everywhere function myidwt(nspat,taut::Array{Float64,4},rhot,t::Array{Float64,4},spatialwlt)
                 wlt = idwt(taut[:,:,1,1] + rhot*t[:,:,1,1],wavelet(spatialwlt[1]))
                     for b in 2:nspat
                         wlt = wlt + idwt(taut[:,:,1,b] + rhot*t[:,:,1,b],wavelet(spatialwlt[b]))
                     end
                 return wlt
            end

@everywhere function mypmap(map)

        taut = map[1]
        t = map[2]
        rhot = map[3]
        x = map[4]
        psf = map[5]
        p = map[6]
        taup = map[7]
        fty = map[8]
        rhop = map[9]
        taus = map[10]
        s = map[11]
        rhos = map[12]
        mu = map[13]
        spatialwlt = map[14]
        μt = map[15]
        nspat = map[16]

        println(myid())
        wlt = myidwt(nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
        b = fty + taup + rhop*p + taus + rhos*s
        wlt_b = wlt + b

        nxy = (size(x))[1]
        nxypsf = (size(psf))[1]
        psfcbe = zeros(Complex64,nxy,nxy)
        psfpad = zeros(Float64,nxy,nxy)
        psfpad[1:nxypsf,1:nxypsf] = psf[:,:]
        psfcbe = 1./(abs(fft(psfpad)).^2+mu)
        x = real(ifft(psfcbe.*fft(wlt_b)))

        # tmp1 = 0.0
        # tmp2 = zeros(Float64,nxy,nxy)
        for b in 1:nspat
                    hx = dwt(x,wavelet(spatialwlt[b]))
                    tmp = hx - taut[:,:,1,b]/rhot
                    t[:,:,1,b] = prox_u(tmp,μt/rhot)
                    taut[:,:,1,b] = taut[:,:,1,b] + rhot*(t[:,:,1,b]-hx)
                    # tmp1 = vecnorm([tmp2 (hx-(t)[:,:,z,b])],2)
                    # tmp2 = (hx-(t)[:,:,z,b])
        end
        # tmp2[:] = 0


        tmp = x-taup/rhop
        p = max(0,tmp)
        taup = taup + rhop*(p-x)

    return wlt,x,t,taut,p,taup
end
nfreq = 1
nspat = 1
rhot = 1.0
μt = 1
mu = 1
rhos = 1
rhop = 1

spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]

@time themap = pmap(mypmap,MAP);
