####################################################################
#######                Main Muffin function                  #######
####################################################################

function muffin_SA(;folder="",dataobj="",datapsf="",nitermax = 500, rhop = 1,
            rhot = 5, rhov = 2, rhos = 1, μt = 5e-1, μv = 1e-0, mueps = 1e-3,
            bw = 25, ws="",parallel="",mask="",dirac="false")

println("")
println("MUFFIN initialisation")

             ##################################
################### data initialisation #################
             ##################################

                 if isempty(dataobj) == false
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
                     elseif dataobj == "halo"
                         psf = "/home/deguignet/simulations/halo/HALO_ONLY_psf.fits"
                         obj = "/home/deguignet/simulations/halo/HALO_ONLY_dirty.fits"
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
                     tmp = string(Pkg.dir("Muffin"))
                     psf = string(tmp,tmp[1],psf)
                     obj = string(tmp,tmp[1],obj)
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
elseif dataobj == "halo"
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
if dirac == "true"
    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,"dirac"]
    println("dirac")
else
    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8]
end

const nspat = length(spatialwlt)
const nfreq = size(psfst.mypsf)[3]
const nspec = 1
const nxy = size(skyst.mydata)[1]
niter = 0
lastiter = 0
#################################

#################################
println("loading param...")
algost = loadparam(nspat,nfreq,nspec,nxy,niter,lastiter,nitermax,spatialwlt)

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
    # admmst = loadarray(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,
    #                     skyst.mydata,psfst.mypsfadj,parallel)
    admmst = loadarray_v2(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,
                        skyst.mydata,psfst.mypsfadj,parallel)
    toolst = loadtools(nitermax,nfreq,nxy)
end

println(typeof(mask))
if typeof(mask) == Int64
    toolst.mask2D = maskgen(admmst.x[:,:,1],mask)
    println(sum(toolst.mask2D)/4)
end

psfst.psfcbe = psfcbe_init(psfst.psfcbe,admmst.x,psfst.mypsf,admmst.mu)

##################################

             ##################################
#####################  Main Admm Loop  #####################
             ##################################
println("Starting ADMM")
#################################
psfst, skyst, algost, admmst, toolst = muffinadmm_SA(psfst, skyst, algost, admmst, toolst)
#################################


return psfst, skyst, algost, admmst, toolst
end



####################################################################
#######                   Main Admm Loop                     #######
####################################################################

function muffinadmm_SA(psfst, skyst, algost, admmst, toolst)

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
const spatialwlt  = algost.spatialwlt
const mask = toolst.mask2D

niter = algost.niter

workers = nworkers()
a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
a = a[1:nfreq]

# println("3D SharedArrays init")
# cube3D = genshared3D(admmst.x,admmst.p,admmst.taup,admmst.wlt,psfst.mypsf,admmst.fty,admmst.s,admmst.taus)
# println("")
# println("4D SharedArrays init")
# cube4D = genshared4D(admmst.t,admmst.taut)

println("3D SharedArrays init")
cube3D = genshared3D_v2(admmst.x,psfst.mypsf,admmst.fty)
println("4D SharedArrays init")
cube4D = genshared4D_v2(nspat,nxy,nfreq)

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


    @sync begin
         @parallel for z in 1:nfreq
              cube3D[z][:,:,4],cube3D[z][:,:,1],cube4D[z][:,:,1,:],cube4D[z][:,:,2,:],
              cube3D[z][:,:,2],cube3D[z][:,:,3] = parallelmuffin(        cube3D[z][:,:,4],
                                                                         cube4D[z][:,:,2,:],
                                                                         cube4D[z][:,:,1,:],rhot,
                                                                         cube3D[z][:,:,1],
                                                                         cube3D[z][:,:,2],
                                                                         cube3D[z][:,:,3],
                                                                         cube3D[z][:,:,5],rhop,
                                                                         cube3D[z][:,:,7],
                                                                         cube3D[z][:,:,6],
                                                                         rhos,mu,spatialwlt,μt,nspat,mask,psfst.psfcbe[:,:,z])
          end
      end

        ##############################
        ######### prox spec ##########
        for z in 1:nfreq
        admmst.s[:,:,z] = copy(cube3D[z][:,:,6])
        admmst.x[:,:,z] = copy(cube3D[z][:,:,1])
        admmst.taus[:,:,z] = copy(cube3D[z][:,:,7])
        admmst.p[:,:,z] = copy(cube3D[z][:,:,2])
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



        for z in 1:nfreq
            cube3D[z][:,:,6] = admmst.s[:,:,z]
            cube3D[z][:,:,7] = admmst.taus[:,:,z]
        end

        ##############################
        ##### computer residues ######

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

#########################################################################

# genshared3D(x,p,taup,wlt,psf,fty,s,taus,listarr)
function genshared3D(x::Array{Float64,3},p::Array{Float64,3},taup::Array{Float64,3},wlt::Array{Float64,3},psf::Array{Float64,3},
                     fty::Array{Float64,3},s::Array{Float64,3},taus::Array{Float64,3})


                    #  d = 2
                    #  x  = x[1:d,1:d,:]
                    #  p =p[1:d,1:d,:]
                    #  taup = taup[1:d,1:d,:]
                    #  wlt =wlt[1:d,1:d,:]
                    #  psf =psf[1:d,1:d,:]
                    #  fty =fty[1:d,1:d,:]
                    #  s =s[1:d,1:d,:]
                    #  taus = taus[1:d,1:d,:]

    Ndim = 7

    nxy = size(x)[1]
    nfreq = size(x)[3]
    workers = nworkers()
    a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
    a = a[1:nfreq]

    tabname = "freq3d"
    listarr = {0 => 1}

    for z in 1:nfreq
        # listarr[z] = SharedArray(Float64,nxy,nxy,Ndim,pids=[1,a[z]])
        listarr[z] = SharedArray(Float64,nxy,nxy,Ndim)

    end

    for z in 1:nfreq
        print("freq","$z"," ")
        listarr[z][:,:,1] = x[:,:,z]
        listarr[z][:,:,2] = p[:,:,z]
        listarr[z][:,:,3] = taup[:,:,z]
        listarr[z][:,:,4] = wlt[:,:,z]
        listarr[z][:,:,5] = fty[:,:,z]
        listarr[z][:,:,6] = s[:,:,z]
        listarr[z][:,:,7] = taus[:,:,z]
    end
    return listarr
end

#########################################################################

# genshared4D(t,taut,listarr)
function genshared4D(t::Array{Float64,4},taut::Array{Float64,4})
    #
    # d = 2
    # t = t[1:d,1:d,:,:]
    # taut = taut[1:d,1:d,:,:]

    listarr = {0 => 1}
    Ndim = 2
    nxy = size(t)[1]
    nfreq = size(t)[3]
    nspat = size(t)[4]
    workers = nworkers()
    a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
    a = a[1:nfreq]

    tabname = "freq4d"
    toeval = ""


    listarr = {0 => 1}

    for z in 1:nfreq
        # listarr[z] = SharedArray(Float64,nxy,nxy,Ndim,nspat,pids=[1,a[z]])
        listarr[z] = SharedArray(Float64,nxy,nxy,Ndim,nspat)
    end

    for z in 1:nfreq
        print("freq","$z"," ")
        for b in 1:nspat
            listarr[z][:,:,1,b] = t[:,:,z,b]
            listarr[z][:,:,2,b] = taut[:,:,z,b]
        end
    end
    return listarr
end

#########################################################################
#########################################################################


function parallelmuffin(wlt::Array{Float64,2},taut::Array{Float64,4},t::Array{Float64,4},rhot::Float64,
                        x::Array{Float64,2},p::Array{Float64,2},taup::Array{Float64,2},
                        fty::Array{Float64,2},rhop::Float64,taus::Array{Float64,2},s::Array{Float64,2},rhos::Float64,
                        mu::Float64,spatialwlt,μt::Float64,nspat::Int,mask::Array{Float64,2},psfcbe::Array{Complex64,2})


    if spatialwlt[end] == "dirac"
        wlt = myidwt_dirac(wlt, nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
        b = fty + taup + rhop*p + taus + rhos*s
        wlt_b = wlt + b

        x = real(ifft(psfcbe.*fft(wlt_b)))

        for b in 1:nspat-1
                    hx = dwt(x,wavelet(spatialwlt[b]))
                    tmp = hx - taut[:,:,1,b]/rhot
                    t[:,:,1,b] = prox_u(tmp,μt/rhot)
                    taut[:,:,1,b] = taut[:,:,1,b] + rhot*(t[:,:,1,b]-hx)
        end
                    hx = x
                    tmp = hx - taut[:,:,1,nspat]/rhot
                    t[:,:,1,nspat] = prox_u(tmp,μt/rhot)
                    taut[:,:,1,nspat] = taut[:,:,1,nspat] + rhot*(t[:,:,1,nspat]-hx)


    else
        wlt = myidwt(wlt, nspat, taut[:,:,1,:], rhot, t[:,:,1,:], spatialwlt)
        b = fty + taup + rhop*p + taus + rhos*s
        wlt_b = wlt + b

        x = real(ifft(psfcbe.*fft(wlt_b)))

        for b in 1:nspat
                    hx = dwt(x,wavelet(spatialwlt[b]))
                    tmp = hx - taut[:,:,1,b]/rhot
                    t[:,:,1,b] = prox_u(tmp,μt/rhot)
                    taut[:,:,1,b] = taut[:,:,1,b] + rhot*(t[:,:,1,b]-hx)
        end
    end



        tmp = x-taup/rhop
        p = max(0,tmp).*mask
        taup = taup + rhop*(p-x)

    return wlt,x,t,taut,p,taup

end


function estime_ssh(s::Array{Float64,3},sh::Array{Float64,3},tmp::Array{Float64,3},
                    nxy::Int64,nspec::Int64,spectralwlt::Array{Float64,3},
                    x::Array{Float64,3},taus::Array{Float64,3},rhov::Float64,rhos::Float64)


    spectralwlt = idct(tmp,3)
    s = (spectralwlt + rhos*x - taus)/(rhov*nspec + rhos)
    sh = dct(s,3)

    return s,sh
end

function loadarray_v2(rhop,rhot,rhov,rhos,μt,μv,mueps,nspat,nfreq,nxy,mydata,mypsfadj,parallel)
    admmst = init_Admmarray(parallel)

    admmst.sh = zeros(Float64,nxy,nxy,nfreq)
    admmst.tauv = zeros(Float64,nxy,nxy,nfreq)
    admmst.v = zeros(Float64,nxy,nxy,nfreq)

        admmst.s = zeros(Float64,nxy,nxy,nfreq)
        admmst.taus = zeros(Float64,nxy,nxy,nfreq)
        admmst.p = zeros(Float64,nxy,nxy,nfreq)
        admmst.x = copy(mydata)
        admmst.fty = cubefilter(mydata,mypsfadj)

    admmst.xmm = zeros(Float64,nxy,nxy,nfreq)
    admmst.spectralwlt = zeros(Float64,nxy,nxy,nfreq)

    admmst.rhop = rhop
    admmst.rhot = rhot
    admmst.rhov = rhov
    admmst.rhos = rhos
    admmst.μt = μt
    admmst.μv = μv
    admmst.mueps = mueps
    admmst.tt = rhot*nspat
    admmst.mu = mueps + rhop + admmst.tt + rhos
    return admmst
end

# cube4D = genshared4D_v2(nspat,nxy,nfreq)
function genshared4D_v2(nspat,nxy,nfreq)

    listarr = {0 => 1}
    Ndim = 2

    workers = nworkers()
    a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
    a = a[1:nfreq]

    listarr = {0 => 1}

    for z in 1:nfreq
        # listarr[z] = SharedArray(Float64,nxy,nxy,Ndim,nspat,pids=[1,a[z]])
        listarr[z] = SharedArray(Float64,nxy,nxy,Ndim,nspat)
    end

    return listarr
end

# cube3D = genshared3D_v2(admmst.x,psfst.mypsf,admmst.fty)
function genshared3D_v2(x::Array{Float64,3},psf::Array{Float64,3},fty::Array{Float64,3})

    Ndim = 7

    nxy = size(x)[1]
    nfreq = size(x)[3]
    workers = nworkers()
    a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
    a = a[1:nfreq]

    listarr = {0 => 1}

    for z in 1:nfreq
        # listarr[z] = SharedArray(Float64,nxy,nxy,Ndim,pids=[1,a[z]])
        listarr[z] = SharedArray(Float64,nxy,nxy,Ndim)

    end

    for z in 1:nfreq
        listarr[z][:,:,1] = x[:,:,z]
        listarr[z][:,:,5] = fty[:,:,z]
    end
    return listarr
end
