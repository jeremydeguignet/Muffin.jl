####################################################################
#######                Main Muffin function                  #######
####################################################################

function muffin_SA(;folder="",dataobj="",datapsf="",nitermax = 500, rhop = 1,
            rhot = 5, rhov = 2, rhos = 1, μt = 5e-1, μv = 1e-0, mueps = 1e-3,
            bw = 1, ws="",parallel="",mask="")

println("")
println("MUFFIN initialisation")

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

println(typeof(mask))
if typeof(mask) == Int64
    toolst.mask2D = maskgen(admmst.x[:,:,1],mask)
    println(sum(toolst.mask2D)/4)
end


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
const mask = toolst.mask2D


spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,WT.haar]

niter = algost.niter

workers = nworkers()
a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
a = a[1:nfreq]

println("4D SharedArrays init")
genshared4D(admmst.t,admmst.taut)
println("3D SharedArrays init")
genshared3D(admmst.x,admmst.p,admmst.taup,admmst.wlt,psfst.mypsf,admmst.fty,admmst.s,admmst.taus)


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
        # for z in 1:nfreq
        #     admmst.wlt[:,:,z],admmst.x[:,:,z],admmst.t[:,:,z,:],admmst.taut[:,:,z,:],admmst.p[:,:,z],admmst.taup[:,:,z] =
        #                                 parallelmuffin(admmst.wlt[:,:,z], admmst.taut[:,:,z,:], admmst.t[:,:,z,:], rhot, admmst.x[:,:,z],
        #                                 psfst.mypsf[:,:,z], admmst.p[:,:,z], admmst.taup[:,:,z],
        #                                 fty[:,:,z], rhop, admmst.taus[:,:,z], admmst.s[:,:,z], rhos, admmst.mu, spatialwlt,
        #                                 μt, nspat,mask)
        #
        # end


            tabname1 = "freq3d"
            tabname2 = "freq4d"

                chaine = string("for z in 1:nfreq;",string(tabname1, "$z", "[:,:,4,:]",",","",
                                                           tabname1, "$z", "[:,:,1,:]",",","",
                                                           tabname2, "$z", "[:,:,1,:]",",","",
                                                           tabname2, "$z", "[:,:,2,:]",",","",
                                                           tabname1, "$z", "[:,:,2,:]",",","",
                                                           tabname1, "$z", "[:,:,3,:]",",","","=",
                                                           "parallelmuffin(",
                                                           tabname1, "$z", "[:,:,4,:]",",","",
                                                           tabname2, "$z", "[:,:,2,:]",",","",
                                                           tabname2, "$z", "[:,:,1,:]",",","","rhot",
                                                           tabname1, "$z", "[:,:,1,:]",",","",
                                                           tabname1, "$z", "[:,:,5,:]",",","",
                                                           tabname1, "$z", "[:,:,2,:]",",","",
                                                           tabname1, "$z", "[:,:,3,:]",",","",
                                                           tabname1, "$z", "[:,:,6,:]",",","","rhop",
                                                           tabname1, "$z", "[:,:,8,:]",",","",
                                                           tabname1, "$z", "[:,:,7,:]",",","",
                                                           "rhos", "mu", "spatialwlt", "μt", "nspat","mask);")," end;")

                toeval = "chaine"

        eval(parse(toeval))


println(eval(parse(toeval)))
        # for z in 1:nfreq
        #
        #     freq3d$z[:,:,4,:],freq3d$z[:,:,1,:],freq4d$z[:,:,1,:],freq4d$z[:,:,2,:],freq3d$z[:,:,2,:],freq3d$z[:,:,3,:] =
        #                                 parallelmuffin(freq3d$z[:,:,4,:],,freq4d$z[:,:,2,:],freq4d$z[:,:,1,:], rhot,freq3d$z[:,:,1,:],
        #                                 freq3d$z[:,:,5,:],freq3d$z[:,:,2,:],freq3d$z[:,:,3,:],freq3d$z[:,:,6,:], rhop,freq3d$z[:,:,8,:],
        #                                 freq3d$z[:,:,7,:], rhos, mu, spatialwlt, μt, nspat,mask)
        #
        #
        #
        # end


        ##############################
        ######### prox spec ##########

        tmp = admmst.tauv + rhov*admmst.v

        admmst.s, admmst.sh = estime_ssh(admmst.s,admmst.sh,tmp,nxy,nspec,admmst.spectralwlt,
                                          admmst.x,admmst.taus,rhov,rhos)

        tmp = admmst.sh - admmst.tauv/rhov
        admmst.v = prox_u(tmp,μv/rhov)


        ########################################
        #### update of Lagrange multipliers ####

        admmst.tauv = admmst.tauv + rhov*(admmst.v-admmst.sh)
        admmst.taus = admmst.taus + rhos*(admmst.s-admmst.x)

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
# Ndim = 8
# listarr = {0 => 1}
# genshared3D(x,p,taup,wlt,psf,fty,s,taus,listarr)

function genshared3D(x::Array{Float64,3},p::Array{Float64,3},taup::Array{Float64,3},wlt::Array{Float64,3},psf::Array{Float64,3},
                     fty::Array{Float64,3},s::Array{Float64,3},taus::Array{Float64,3})


    Ndim = 8
    listarr = {0 => 1}
    nxy = size(x)[1]
    nfreq = size(x)[3]
    workers = nworkers()
    a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
    a = a[1:nfreq]

    tabname = "freq3d"
    toeval = ""

    for z in 1:nfreq
        toeval = string(toeval,string(tabname, "$z", "=", "SharedArray(Float64,$nxy,$nxy,$Ndim,pids=[1,$a[$z]]);"))
    end
    eval(parse(toeval))

    toeval = ""
    result = {0 => 1}
    for z in 1:nfreq
        result[1] = x[:,:,z]
        result[2] = p[:,:,z]
        result[3] = taup[:,:,z]
        result[4] = wlt[:,:,z]
        result[5] = psf[:,:,z]
        result[6] = fty[:,:,z]
        result[7] = s[:,:,z]
        result[8] = taus[:,:,z]
        chaine = string("for nd in 1:$Ndim;",string(tabname, "$z", "[:,:,nd]","=","$result[nd];")," end;")
        toeval = string(toeval,chaine)
    end
    eval(parse(toeval))
end


#########################################################################
# Ndim = 2
# listarr = {0 => 1}
# genshared4D(t,taut,listarr)

function genshared4D(t::Array{Float64,4},taut::Array{Float64,4})

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

    for z in 1:nfreq
        toeval = string(toeval,string(tabname, "$z", "=", "SharedArray(Float64,$nxy,$nxy,$Ndim,$nspat,pids=[1,$a[$z]]);"))
    end
    eval(parse(toeval))

    toeval = ""
    result = {0 => 1}
    for z in 1:nfreq
        for b in 1:nspat
            result[1] = t[:,:,z,b]
            result[2] = taut[:,:,z,b]
            chaine = string("for nd in 1:$Ndim;",string(tabname, "$z", "[:,:,nd,$b]","=","$result[nd];")," end;")
            toeval = string(toeval,chaine)
        end
    end
    eval(parse(toeval))
end

#########################################################################
# nxy = 2
# nfreq = 4
# nspat = 8
#
# workers = nworkers()
# a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
# a = a[1:nfreq]
#
# x = randn(nxy,nxy,nfreq);
# p =randn(nxy,nxy,nfreq);
# taup = randn(nxy,nxy,nfreq);
# wlt = randn(nxy,nxy,nfreq);
# psf =randn(nxy,nxy,nfreq);
# fty = randn(nxy,nxy,nfreq);
# s =randn(nxy,nxy,nfreq);
# taus =randn(nxy,nxy,nfreq);
# t =randn(nxy,nxy,nfreq,nspat);
# taut =randn(nxy,nxy,nfreq,nspat);

# #########################################################################
# tabname = "freq"
# toeval = ""
#
# listarr = {0 => 1}
# for z in 1:nfreq
#       listarr[z] = SharedArray(Float64,nxy,nxy,Ndim,pids=[1,a[z]]);
# end
#
# for n in 1:nfreq
#     proc = n
#   toeval = string(toeval,string(tabname, "$n", "=", "similar(listarr[$n]);"))
# end
# eval(parse(toeval))
# toeval = ""
#
# result = {0 => 1}
# for z in 1:nfreq
#
#     result[1] = x[:,:,z]
#     result[2] = p[:,:,z]
#     result[3] = taup[:,:,z]
#     result[4] = wlt[:,:,z]
#     result[5] = psf[:,:,z]
#     result[6] = fty[:,:,z]
#     result[7] = s[:,:,z]
#     result[8] = taus[:,:,z]
#
#     chaine = string("for nd in 1:Ndim;",string(tabname, "$z", "[:,:,nd]","=","$result[nd];")," end;")
#
#     toeval = string(toeval,chaine)
# end
# eval(parse(toeval))
#########################################################################
