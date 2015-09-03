# dim = 1000000000000000000000000
#
# a = @spawnat(2,ones(dim,dim))
# b = @spawnat(2,2*ones(dim,dim))
# @time c = @spawnat(2,fetch(a)+fetch(b))
#
# image = randn(2048,2048,4);
# psf = randn(2048,2048,4);
#
# imlist = [im1, im2, im3, im4]
# [imlist[z] = @spawnat(z+1,image[:,:,z]) for z in 1:4]
#
# @everywhere using Images
#
# im1 = @spawnat(2,image[:,:,1])
# im2 = @spawnat(3,image[:,:,2])
# im3 = @spawnat(4,image[:,:,3])
# im4 = @spawnat(5,image[:,:,4])
#
# psf1 = @spawnat(2,psf[:,:,1])
# psf2 = @spawnat(3,psf[:,:,2])
# psf3 = @spawnat(4,psf[:,:,3])
# psf4 = @spawnat(5,psf[:,:,4])
#
# resultat1 = @spawnat(2,conv2(fetch(im1),fetch(psf1)))
# resultat2 = @spawnat(3,conv2(fetch(im2),fetch(psf2)))
# resultat3 = @spawnat(4,conv2(fetch(im3),fetch(psf3)))
# resultat4 = @spawnat(5,conv2(fetch(im4),fetch(psf4)))
#
#
# n = size(image)[3]
# psf = Array(ASCIIString,n)
# img = Array(ASCIIString,n)
# resultat =  Array(ASCIIString,n)
# for z in 1:n
#     psf[z] = string("psf",z)
#     img[z] = string("img",z)
#     resultat[z] = string("resultat",z)
# end
#
# for z in 1:n
#     ex = parse ("  toto = 2x + y ^ 2 + 1")
#     psf$z = (eval(ex))
# end

####################################################
function cubetoproc(cube::Array{Float64,3})
    nfreq = size(cube)[3]
    workers = nworkers()


    if workers >= nfreq
        cuberef = Array(RemoteRef,nfreq)
        for z in 1:nfreq
        cuberef[z] = @spawnat(z+1,cube[:,:,z])
        end
    else
        a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
        a = a[1:nfreq]
        cuberef = Array(RemoteRef,nfreq)
        for z in 1:nfreq
            cuberef[z] = @spawnat(a[z],cube[:,:,z])
        end

    end

    return cuberef

end

function cubetoproc(cube::Array{Float64,4})
    nfreq = size(cube)[3]
    workers = nworkers()


    if workers >= nfreq
        cuberef = Array(RemoteRef,nfreq)
        for z in 1:nfreq
        cuberef[z] = @spawnat(z+1,cube[:,:,z,:])
        end
    else
        a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
        a = a[1:nfreq]
        cuberef = Array(RemoteRef,nfreq)
        for z in 1:nfreq
        cuberef[z] = @spawnat(a[z],cube[:,:,z,:])
        end

    end

    return cuberef

end

type REF
    wlt::Array{RemoteRef}
    taut::Array{RemoteRef}
    t::Array{RemoteRef}
    x::Array{RemoteRef}
    psf::Array{RemoteRef}
    p::Array{RemoteRef}
    taup::Array{RemoteRef}
    fty::Array{RemoteRef}
    taus::Array{RemoteRef}
    s::Array{RemoteRef}
end

function init_REF()
    return REF([],[],[],[],[],[],[],[],[],[])
end

function loadref(wlt::Array{Float64,3}, taut::Array{Float64,4}, t::Array{Float64,4}, x::Array{Float64,3}, psf::Array{Float64,3},
                     p::Array{Float64,3}, taup::Array{Float64,3}, fty::Array{Float64,3}, taus::Array{Float64,3},s::Array{Float64,3})
    refst = init_REF()
    refst.wlt = cubetoproc(wlt)
    refst.taut = cubetoproc(taut)
    refst.t = cubetoproc(t)
    refst.x = cubetoproc(x)
    refst.psf = cubetoproc(psf)
    refst.p = cubetoproc(p)
    refst.taup = cubetoproc(taup)
    refst.fty = cubetoproc(fty)
    refst.taus = cubetoproc(taus)
    refst.s = cubetoproc(s)

    return refst

end
#
# image = randn(2048,2048,4);
# psf = randn(2048,2048,4);
#
# im = cubetoproc(image);
# dirt = cubetoproc(psf);



function muffin_DA(;folder="",dataobj="",datapsf="",nitermax = 500, rhop = 1,
                rhot = 5, rhov = 2, rhos = 1, μt = 5e-1, μv = 1e-0, mueps = 1e-3,
                bw = 5, ws="",parallel="")

println("")
println("MUFFIN_DA initialisation")

                 ##################################
    ################### data initialisation #################
                 ##################################


    if typeof(dataobj) == ASCIIString
        if dataobj == "m31"
            psf = "data/meerkat_m30_25pix.psf.fits"
            obj = "data/M31.fits"
            tmp = string(Pkg.dir("Muffin"))
            psf = string(tmp,tmp[1],psf)
            obj = string(tmp,tmp[1],obj)
        elseif dataobj == "andro"
            psf = "data/meerkat_m30_25pix.psf.fits"
            obj = "data/andro.fits"
            tmp = string(Pkg.dir("Muffin"))
            psf = string(tmp,tmp[1],psf)
            obj = string(tmp,tmp[1],obj)
        elseif dataobj == "2gauss"
            psf = "data/meerkat_m30_25pix.psf.fits"
            obj = "data/2gauss.fits"
            tmp = string(Pkg.dir("Muffin"))
            psf = string(tmp,tmp[1],psf)
            obj = string(tmp,tmp[1],obj)
        elseif dataobj == "chiara"
            psf = "/home/deguignet/Julia/example_sim_psf.fits"
            obj = "/home/deguignet/Julia/example_sim_dirty.fits"
        elseif isempty(folder)
            tmp = pwd()
            psf = string(tmp,tmp[1],datapsf)
            obj = string(tmp,tmp[1],dataobj)
        elseif typeof(folder) == ASCIIString
                 psf = string(folder,folder[1],datapsf)
                 obj = string(folder,folder[1],dataobj)
        else
                 error("data folder is not correct")
        end

    elseif isempty(dataobj)
        psf = "data/meerkat_m30_25pix.psf.fits"
        obj = "data/M31.fits"
    end

println("psf :"," ",psf)
println("obj :"," ",obj)

                ##################################
    ################# Structure initialisation #################
                ##################################

if dataobj == "chiara"
    ##################################
    println("loading psf...")
    psfst = loadpsf_dirty(psf)
    println("loading sky...")
    skyst = loadsky_dirty(obj,psfst.mypsf,psfst.nu)
    ##################################
else
    ##################################
    println("loading psf...")
    psfst = loadpsf(psf,bw)
    println("loading sky...")
    skyst = loadsky(obj,psfst.mypsf,psfst.nu)
    ##################################
end


    ##################################
    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8]
    # spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]

    const nspat = length(spatialwlt)
    const nfreq = size(psfst.mypsf)[3]
    const nspec = 1
    const nxy = size(skyst.mydata)[1]
    niter = 0
    lastiter = 0
    #################################

    #################################
    println("loading param...")
    algost = loadparam(nspat,nfreq,nspec,nxy,niter,lastiter,nitermax)

    println("loading arrays...")
    if ws == "true"
        file = string("/home/deguignet/resultats_400ite_woshar.jld")
        # file = string("/Users/deguignet/test.jld")
        admmst = loadarray_ws(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,
                            skyst.mydata,psfst.mypsfadj,file)
    algost.niter = load(file,"algost.lastiter")
    toolst = loadtools(nitermax,nfreq,nxy)
    toolst.tol1 = load(file,"toolst.tol1")
    toolst.tol2 = load(file,"toolst.tol2")
    toolst.tol3 = load(file,"toolst.tol3")
    toolst.tol4 = load(file,"toolst.tol4")
    toolst.tol5 = load(file,"toolst.tol5")
    skyst.mydata = load(file,"skyst.mydata")


    else
        admmst = loadarray(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,
                            skyst.mydata,psfst.mypsfadj,parallel)
        toolst = loadtools(nitermax,nfreq,nxy)
    end

    println("loading ref...")
    refst = loadref(admmst.wlt, admmst.taut, admmst.t, admmst.x, psfst.mypsf, admmst.p, admmst.taup, admmst.fty, admmst.taus, admmst.s)

    ##################################

                 ##################################
    #####################  Main Admm Loop  #####################
                 ##################################
    println("Starting ADMM")
    #################################
    psfst, skyst, algost, admmst, toolst = muffinadmm_DA(psfst, skyst, algost, admmst, toolst, refst)
    #################################


    return psfst, skyst, algost, admmst, toolst
end



function muffinadmm_DA(psfst, skyst, algost, admmst, toolst, refst)

    const rhop = admmst.rhop
    const rhot = admmst.rhot
    const rhov = admmst.rhov
    const rhos = admmst.rhos
    const μt = admmst.μt
    const μv = admmst.μv
    const mueps = admmst.mueps
    const tt = admmst.tt
    const mu = admmst.mu
    const nspat = algost.nspat
    const nfreq = algost.nfreq
    const nspec = algost.nspec
    const nxy = algost.nxy
    const fty = admmst.fty
    const nitermax = algost.nitermax

    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]

    niter = algost.niter

    workers = nworkers()
    a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
    a = a[1:nfreq]

    loop = true

    println("")
    println("                                =======================                                ")
    println("                               | MUFFIN RECONSTRUCTION |                               ")

    println("========================================================================================")
    println("ADMM"," ","|"," ","||x - xm||"," ","|"," ","||x - xp||"," ","|"," ",
            "||Hx - t||"," ","|"," ","||x - s||"," "," |"," ","||sh - v||"," ","|"," ","     time      |")
    println("========================================================================================")
    tic()
        while loop
            tic()
            niter +=1
            println("    "," ","|"," ","          "," ","|"," ","          "," ","|"," ",
                    "          "," ","|"," ","         "," "," |"," ","          "," ","|"," ")
            # println("ADMM iteration: $niter")

    ######################################
    ######################################


    @sync @parallel for z in 1:nfreq
        # # refst.wlt[z],refst.x[z],refst.t[z],refst.taut[z],refst.p[z],refst.taup[z] =
        #             println(                @spawnat(a[z],parallelmuffin(fetch(refst.wlt[z]), fetch(refst.taut[z]), fetch(refst.t[z]), rhot, fetch(refst.x[z]),
        #                             fetch(refst.psf[z]), fetch(refst.p[z]), fetch(refst.taup[z]),
        #                             fetch(refst.fty[z]), rhop, fetch(refst.taus[z]), fetch(refst.s[z]), rhos, admmst.mu, spatialwlt,
        #                             μt, nspat)))
        refst.wlt[z] =  @spawnat(a[z],muffpar_wlt(fetch(refst.wlt[z]), fetch(refst.taut[z]), fetch(refst.t[z]), rhot,  spatialwlt, nspat))
        refst.x[z] = @spawnat(a[z],muffpar_x(fetch(refst.wlt[z]), fetch(refst.x[z]), fetch(refst.psf[z]), fetch(refst.p[z]),
                                             fetch(refst.taup[z]), fetch(refst.fty[z]), rhop, admmst.taus[:,:,z], admmst.s[:,:,z],
                                             rhos, admmst.mu, nspat))
        refst.t[z] = @spawnat(a[z],muffpar_t(fetch(refst.taut[z]), fetch(refst.t[z]), rhot, fetch(refst.x[z]), spatialwlt, μt, nspat))
        refst.taut[z] = @spawnat(a[z],muffpar_taut(fetch(refst.taut[z]), fetch(refst.t[z]), rhot, fetch(refst.x[z]), spatialwlt, μt, nspat))
        refst.p[z] = @spawnat(a[z],muffpar_p(fetch(refst.x[z]), fetch(refst.p[z]), fetch(refst.taup[z]), rhop))
        refst.taup[z] = @spawnat(a[z],muffpar_taup(fetch(refst.x[z]), fetch(refst.p[z]), fetch(refst.taup[z]), rhop))


    end

    for z in 1: nfreq
        admmst.x[:,:,z] = fetch(refst.x[z])
        admmst.p[:,:,z] = fetch(refst.p[z])
    end

    tmp = admmst.tauv + rhov*admmst.v

    admmst.s, admmst.sh = estime_ssh(admmst.s,admmst.sh,tmp,nxy,nspec,admmst.spectralwlt,
                                      admmst.x,admmst.taus,rhov,rhos)

    tmp = admmst.sh - admmst.tauv/rhov
    admmst.v = prox_u(tmp,μv/rhov)


    ########################################
    #### update of Lagrange multipliers ####

    admmst.tauv = admmst.tauv + rhov*(admmst.v-admmst.sh)
    admmst.taus = admmst.taus + rhos*(admmst.s-admmst.x)

    ######################################
    ######################################

            ##############################
            ##### computer residues ######

            # push!(toolst.tol1,vecnorm(admmst.x - admmst.xmm, 2)^2)
            # push!(toolst.tol2,vecnorm(admmst.x - admmst.p, 2)^2)
            # # push!(toolst.tol3,tmp1^2)
            # push!(toolst.tol4,vecnorm(admmst.x - admmst.s, 2)^2)
            # push!(toolst.tol5,vecnorm(admmst.sh - admmst.v, 2)^2)

            push!(toolst.tol1,vecnorm(admmst.x - admmst.xmm, 2)^2)
            push!(toolst.tol2,vecnorm(admmst.x - admmst.p, 2)^2)
            # push!(toolst.tol3,tmp1^2)
            push!(toolst.tol4,vecnorm(admmst.x - admmst.s, 2)^2)
            push!(toolst.tol5,vecnorm(admmst.sh - admmst.v, 2)^2)


            # ##############################
            # ####### stopping rule ########

            if (niter >= nitermax) #|| ((toolst.tol1[niter] < 1E-6) && (toolst.tol2[niter] < 1E-4))
                loop = false
                algost.lastiter = niter
            end

            admmst.xmm[:] = admmst.x

            @printf("%03d  | ", niter)
            @printf("%02.04e | ", toolst.tol1[niter])
            @printf("%02.04e | ", toolst.tol2[niter])
            # @printf("%02.04e | ", toolst.tol3[niter])
            @printf("%02.04e | ", 0)
            @printf("%02.04e | ", toolst.tol4[niter])
            @printf("%02.04e | ", toolst.tol5[niter])
            @printf("%f seconds  \n", toq())

        end
    println("")
    @printf("time for ADMM : %f seconds \n", toq())
    ####################################################################
    ####################################################################
    ####################################################################

    return psfst, skyst, algost, admmst, toolst

end








function muffpar_wlt(wlt::Array{Float64,2},taut::Array{Float64,4},t::Array{Float64,4},rhot::Float64,
                     spatialwlt, nspat::Int)
    wlt = myidwt(wlt, nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
end
function muffpar_x(wlt::Array{Float64,2}, x::Array{Float64,2},psf::Array{Float64,2},p::Array{Float64,2},taup::Array{Float64,2},
                        fty::Array{Float64,2},rhop::Float64,taus::Array{Float64,2},s::Array{Float64,2},rhos::Float64,
                        mu::Float64,nspat::Int)
    b = fty + taup + rhop*p + taus + rhos*s
    wlt_b = wlt + b

    nxy = (size(x))[1]
    nxypsf = (size(psf))[1]
    psfcbe = zeros(Complex64,nxy,nxy)
    psfpad = zeros(Float64,nxy,nxy)
    psfpad[1:nxypsf,1:nxypsf] = psf[:,:]
    psfcbe = 1./(abs(fft(psfpad)).^2+mu)
    x = real(ifft(psfcbe.*fft(wlt_b)))
end
function muffpar_t(taut::Array{Float64,4},t::Array{Float64,4},rhot::Float64,
                        x::Array{Float64,2},spatialwlt,μt::Float64,nspat::Int)
    for b in 1:nspat
                hx = dwt(x,wavelet(spatialwlt[b]))
                tmp = hx - taut[:,:,1,b]/rhot
                t[:,:,1,b] = prox_u(tmp,μt/rhot)
    end
    return t
end
function muffpar_taut(taut::Array{Float64,4},t::Array{Float64,4},rhot::Float64,
                        x::Array{Float64,2},spatialwlt,μt::Float64,nspat::Int)
    for b in 1:nspat
                hx = dwt(x,wavelet(spatialwlt[b]))
                tmp = hx - taut[:,:,1,b]/rhot
                t[:,:,1,b] = prox_u(tmp,μt/rhot)
                taut[:,:,1,b] = taut[:,:,1,b] + rhot*(t[:,:,1,b]-hx)
    end
    return taut
end
function muffpar_p(x::Array{Float64,2},p::Array{Float64,2},taup::Array{Float64,2},rhop::Float64)

            tmp = x-taup/rhop
            p = max(0,tmp)
end
function muffpar_taup(x::Array{Float64,2},p::Array{Float64,2},taup::Array{Float64,2},rhop::Float64)
        taup = taup + rhop*(p-x)
end
