####################################################################
#######                Main Muffin function                  #######
####################################################################

function muffin_random(;folder="",dataobj="",datapsf="",nitermax = 500, rhop = 1,
            rhot = 5, rhov = 2, rhos = 1, μt = 5e-1, μv = 1e-0, mueps = 1e-3,
            bw = 25, ws="",parallel="",mask="",dirac="false")

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
if dirac == "true"
    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8,"dirac"]
    println("dirac")
else
    spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8]
end
# spatialwlt  = [WT.db1,WT.db2,WT.db3,WT.db4,WT.db5,WT.db6,WT.db7,WT.db8]
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
psfst, skyst, algost, admmst, toolst = muffinadmm_random(psfst, skyst, algost, admmst, toolst)
#################################


return psfst, skyst, algost, admmst, toolst
end



####################################################################
#######                   Main Admm Loop                     #######
####################################################################

function muffinadmm_random(psfst, skyst, algost, admmst, toolst)

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
# const mask = toolst.mask2D[1:2,1:2,:]





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
println("")
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
                ###########################
                ### Ranndom permutation ###
                ###########################

                permutation = randperm(5)
        
                for p in 1:5
                    perm = p

                    if perm == 1
                        @sync begin
                             @parallel for z in 1:nfreq
                                  cube3D[z][:,:,1] = permuteX( cube3D[z][:,:,4],
                                                               cube4D[z][:,:,2,:],
                                                               cube4D[z][:,:,1,:],rhot,
                                                               cube3D[z][:,:,1],
                                                               cube3D[z][:,:,2],
                                                               cube3D[z][:,:,3],
                                                               cube3D[z][:,:,6],rhop,
                                                               cube3D[z][:,:,8],
                                                               cube3D[z][:,:,7],
                                                               rhos,spatialwlt,nspat,psfst.psfcbe[:,:,z])
                              end
                          end

                    elseif perm == 2
                        @sync begin
                             @parallel for z in 1:nfreq
                                  cube4D[z][:,:,1,:] = permuteT(cube4D[z][:,:,2,:],cube4D[z][:,:,1,:],rhot,cube3D[z][:,:,1],spatialwlt,μt,nspat)
                              end
                          end

                    elseif perm == 3
                        @sync begin
                             @parallel for z in 1:nfreq
                                  cube3D[z][:,:,2] = permuteP( cube3D[z][:,:,1], cube3D[z][:,:,2], cube3D[z][:,:,3], rhop, mask)
                              end
                          end

                    elseif perm == 4

                        for z in 1:nfreq
                        admmst.s[:,:,z] = copy(cube3D[z][:,:,7])
                        admmst.x[:,:,z] = copy(cube3D[z][:,:,1])
                        admmst.taus[:,:,z] = copy(cube3D[z][:,:,8])
                        end

                        tmp = admmst.tauv + rhov*admmst.v

                        admmst.s, admmst.sh = estime_ssh(admmst.s,admmst.sh,tmp,nxy,nspec,admmst.spectralwlt,
                                                          admmst.x,admmst.taus,rhov,rhos)

                        for z in 1:nfreq
                            cube3D[z][:,:,7] = admmst.s[:,:,z]
                        end

                    elseif perm == 5
                        tmp = admmst.sh - admmst.tauv/rhov
                        admmst.v = prox_u(tmp,μv/rhov)
                    end
                end




        ########################################
        #### update of Lagrange multipliers ####


        admmst.tauv = admmst.tauv + rhov*(admmst.v-admmst.sh)
        admmst.taus = admmst.taus + rhos*(admmst.s-admmst.x)



        for z in 1:nfreq
            cube3D[z][:,:,3] = cube3D[z][:,:,3] + rhop*(cube3D[z][:,:,2]-cube3D[z][:,:,1]) # maj taup
            cube3D[z][:,:,8] = admmst.taus[:,:,z]
            admmst.p[:,:,z] = copy(cube3D[z][:,:,2])
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
