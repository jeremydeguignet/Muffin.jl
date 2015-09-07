####################################################################
#######                Main Muffin function                  #######
####################################################################

function muffin_SA(;folder="",dataobj="",datapsf="",nitermax = 500, rhop = 1,
            rhot = 5, rhov = 2, rhos = 1, μt = 5e-1, μv = 1e-0, mueps = 1e-3,
            bw = 5, ws="",parallel="",mask="")

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
        # tic()
        # @parallel for i in 2:nfreq+1
        #     z = i-1
        #     admmst.wlt[:,:,z],admmst.x[:,:,z],admmst.t[:,:,z,:],admmst.taut[:,:,z,:],admmst.p[:,:,z],admmst.taup[:,:,z] =
        #                                 @fetchfrom(i,parallelmuffin(admmst.wlt[:,:,z], admmst.taut[:,:,z,:], admmst.t[:,:,z,:], rhot, admmst.x[:,:,z],
        #                                 psfst.mypsf[:,:,z], admmst.p[:,:,z], admmst.taup[:,:,z],
        #                                 fty[:,:,z], rhop, admmst.taus[:,:,z], admmst.s[:,:,z], rhos, admmst.mu, spatialwlt,
        #                                 μt, nspat))
########################

        for z in 1:nfreq
            admmst.wlt[:,:,z],admmst.x[:,:,z],admmst.t[:,:,z,:],admmst.taut[:,:,z,:],admmst.p[:,:,z],admmst.taup[:,:,z] =
                                        parallelmuffin(admmst.wlt[:,:,z], admmst.taut[:,:,z,:], admmst.t[:,:,z,:], rhot, admmst.x[:,:,z],
                                        psfst.mypsf[:,:,z], admmst.p[:,:,z], admmst.taup[:,:,z],
                                        fty[:,:,z], rhop, admmst.taus[:,:,z], admmst.s[:,:,z], rhos, admmst.mu, spatialwlt,
                                        μt, nspat,mask)

        end






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

######################################
######################################

        #
        # ##############################
        # ########## update x ##########
        # for z in 1:nfreq
        #     admmst.wlt[:,:,z] = myidwt((admmst.wlt)[:,:,z], nspat, (admmst.taut)[:,:,z,:], rhot,
        #                         (admmst.t)[:,:,z,:], spatialwlt)
        # end
        #
        # b = fty + admmst.taup + rhop*admmst.p + admmst.taus + rhos*admmst.s
        #
        # admmst.x = estime_x_par(admmst.x,psfst.mypsf,admmst.wlt + b,mu,nfreq)
        # ##############################
        # ######### prox spat ##########
        # tmp1 = 0.0
        # tmp2 = zeros(Float64,nxy,nxy)
        # for z in 1:nfreq
        #     for b in 1:nspat
        #             hx = dwt(admmst.x[:,:,z],wavelet(spatialwlt[b]))
        #             tmp = hx - admmst.taut[:,:,z,b]/rhot
        #             admmst.t[:,:,z,b] = prox_u(tmp,μt/rhot)
        #             admmst.taut[:,:,z,b] = admmst.taut[:,:,z,b] + rhot*(admmst.t[:,:,z,b]-hx)
        #             tmp1 = vecnorm([tmp2 (hx-(admmst.t)[:,:,z,b])],2)
        #             tmp2 = (hx-(admmst.t)[:,:,z,b])
        #     end
        # end
        # tmp2[:] = 0
        # ##############################
        # ###### prox positivity #######
        #
        # tmp = admmst.x-admmst.taup/rhop
        #
        # admmst.p = max(0,tmp)
        #
        # ##############################
        # ######### prox spec ##########
        #
        # tmp = admmst.tauv + rhov*admmst.v
        #
        # admmst.s, admmst.sh = estime_ssh(admmst.s,admmst.sh,tmp,nxy,nspec,admmst.spectralwlt,
        #                                   admmst.x,admmst.taus,rhov,rhos)
        #
        # tmp = admmst.sh - admmst.tauv/rhov
        # admmst.v = prox_u(tmp,μv/rhov)
        #
        #
        # ########################################
        # #### update of Lagrange multipliers ####
        #
        # admmst.taup = admmst.taup + rhop*(admmst.p-admmst.x)
        # admmst.tauv = admmst.tauv + rhov*(admmst.v-admmst.sh)
        # admmst.taus = admmst.taus + rhos*(admmst.s-admmst.x)

######################################
######################################


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

        #
        # println("ADMM"," ","|"," ","||x - xm||"," ","|"," ","||x - xp||"," ","|"," ",
        #         "||Hx - t||"," ","|"," ","||x - s||"," ","|"," ","||sh - v||"," ","|"," ","time")
        # println(niter," ","|"," ",toolst.tol1[niter]," ","|"," ",toolst.tol2[niter]," ","|"," ",
        #         toolst.tol3[niter]," ","|"," ",toolst.tol4[niter]," ","|"," ",toolst.tol5[niter]," ","|"," ",toq())



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

function cubetosharedA(x::Array{Float64,3},p::Array{Float64,3},taup::Array{Float64,3},wlt::Array{Float64,3},psf::Array{Float64,3},
                       fty::Array{Float64,3},s::Array{Float64,3},taus::Array{Float64,3})

    nfreq = size(x)[3]
    nxy = size(x)[1]

    workers = nworkers()
    a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
    a = a[1:nfreq]

    # -------------------------------------

    tabname = "SA_freq"
    toeval = ""
    result = Array(Float64,nxy,nxy,8)
    for z in 1:nfreq
        for b in 1:8

        result[:,:,1] = x[:,:,z]
        result[:,:,2] = p[:,:,z]
        result[:,:,3] = taup[:,:,z]
        result[:,:,4] = wlt[:,:,z]
        result[:,:,5] = psf[:,:,z]
        result[:,:,6] = fty[:,:,z]
        result[:,:,7] = s[:,:,z]
        result[:,:,8] = taus[:,:,z]
        toeval = string(toeval,string(tabname,"$z","_","$b","=","result[:,:,$b];"))
        end
    end
    eval(parse(toeval))

    for z in 1:nfreq


    for z in 1:nfreq
        SA_freq[z] = SA_freq$z
    end



    return toeval


end

d = 2

x = randn(d,d,4);
p =randn(d,d,4);
taup = randn(d,d,4);
wlt = randn(d,d,4);
psf =randn(d,d,4);
fty = randn(d,d,4);
s =randn(d,d,4);
taus =randn(d,d,4);

# result = SharedArray(Float64,nxy,nxy,8,pids=[1,a[z]])
result = 1
# result[:,:,1] = x[:,:,z]
# result[:,:,2] = p[:,:,z]
# result[:,:,3] = taup[:,:,z]
# result[:,:,4] = wlt[:,:,z]
# result[:,:,5] = psf[:,:,z]
# result[:,:,6] = fty[:,:,z]
# result[:,:,7] = s[:,:,z]
# result[:,:,8] = taus[:,:,z]

for z in 1:nfreq
    result = cat(3,x[:,:,z],p[:,:,z],taup[:,:,z],wlt[:,:,z],psf[:,:,z],fty[:,:,z],s[:,:,z],taus[:,:,z])
    toeval = string(toeval,string(tabname,"$z","=","$result;"))
end

function genshared3D(x::Array{Float64,2},p::Array{Float64,2},taup::Array{Float64,2},wlt::Array{Float64,2},psf::Array{Float64,2},
                     fty::Array{Float64,2},s::Array{Float64,2},taus::Array{Float64,2},proc::Int64)

    nxy = size(x)[1]

    cube = SharedArray(Float64,nxy,nxy,8,pids=[1,proc])
    cube[:,:,1] = x
    cube[:,:,2] = p
    cube[:,:,3] = taup
    cube[:,:,4] = wlt
    cube[:,:,5] = psf
    cube[:,:,6] = fty
    cube[:,:,7] = s
    cube[:,:,8] = taus

    return cube
end

function genshared4D(t::Array{Float64,4},taut::Array{Float64,4},nspat::Int64,proc::Int64)

    nxy = size(x)[1]

    cube = SharedArray(Float64,nxy,nxy,2,nspat,8,pids=[1,proc])
    cube[:,:,1,:] = t
    cube[:,:,2,:] = taut

    return cube
end

freq1 = genshared3D(admmst.x[:,:,1],admmst.p[:,:,1],admmst.taup[:,:,1],admmst.wlt[:,:,1],psfst.mypsf[:,:,1],admmst.fty[:,:,1],admmst.s[:,:,1],admmst.taus[:,:,1],a[1])
freq2 = genshared3D(admmst.x[:,:,2],admmst.p[:,:,2],admmst.taup[:,:,2],admmst.wlt[:,:,2],psfst.mypsf[:,:,2],admmst.fty[:,:,2],admmst.s[:,:,2],admmst.taus[:,:,2],a[2])
freq3 = genshared3D(admmst.x[:,:,3],admmst.p[:,:,3],admmst.taup[:,:,3],admmst.wlt[:,:,3],psfst.mypsf[:,:,3],admmst.fty[:,:,3],admmst.s[:,:,3],admmst.taus[:,:,3],a[3])
freq4 = genshared3D(admmst.x[:,:,4],admmst.p[:,:,4],admmst.taup[:,:,4],admmst.wlt[:,:,4],psfst.mypsf[:,:,4],admmst.fty[:,:,4],admmst.s[:,:,4],admmst.taus[:,:,4],a[4])
freq5 = genshared3D(admmst.x[:,:,5],admmst.p[:,:,5],admmst.taup[:,:,5],admmst.wlt[:,:,5],psfst.mypsf[:,:,5],admmst.fty[:,:,5],admmst.s[:,:,5],admmst.taus[:,:,5],a[5])
freq6 = genshared3D(admmst.x[:,:,6],admmst.p[:,:,6],admmst.taup[:,:,6],admmst.wlt[:,:,6],psfst.mypsf[:,:,6],admmst.fty[:,:,6],admmst.s[:,:,6],admmst.taus[:,:,6],a[6])
freq7 = genshared3D(admmst.x[:,:,7],admmst.p[:,:,7],admmst.taup[:,:,7],admmst.wlt[:,:,7],psfst.mypsf[:,:,7],admmst.fty[:,:,7],admmst.s[:,:,7],admmst.taus[:,:,7],a[7])
freq8 = genshared3D(admmst.x[:,:,8],admmst.p[:,:,8],admmst.taup[:,:,8],admmst.wlt[:,:,8],psfst.mypsf[:,:,8],admmst.fty[:,:,8],admmst.s[:,:,8],admmst.taus[:,:,8],a[8])
freq9 = genshared3D(admmst.x[:,:,9],admmst.p[:,:,9],admmst.taup[:,:,9],admmst.wlt[:,:,9],psfst.mypsf[:,:,9],admmst.fty[:,:,9],admmst.s[:,:,9],admmst.taus[:,:,9],a[9])
freq10 = genshared3D(admmst.x[:,:,10],admmst.p[:,:,10],admmst.taup[:,:,10],admmst.wlt[:,:,10],psfst.mypsf[:,:,10],admmst.fty[:,:,10],admmst.s[:,:,10],admmst.taus[:,:,10],a[10])
freq11 = genshared3D(admmst.x[:,:,11],admmst.p[:,:,11],admmst.taup[:,:,11],admmst.wlt[:,:,11],psfst.mypsf[:,:,11],admmst.fty[:,:,11],admmst.s[:,:,11],admmst.taus[:,:,11],a[11])
freq12 = genshared3D(admmst.x[:,:,12],admmst.p[:,:,12],admmst.taup[:,:,12],admmst.wlt[:,:,12],psfst.mypsf[:,:,12],admmst.fty[:,:,12],admmst.s[:,:,12],admmst.taus[:,:,12],a[12])
freq13 = genshared3D(admmst.x[:,:,13],admmst.p[:,:,13],admmst.taup[:,:,13],admmst.wlt[:,:,13],psfst.mypsf[:,:,13],admmst.fty[:,:,13],admmst.s[:,:,13],admmst.taus[:,:,13],a[13])
freq14 = genshared3D(admmst.x[:,:,14],admmst.p[:,:,14],admmst.taup[:,:,14],admmst.wlt[:,:,14],psfst.mypsf[:,:,14],admmst.fty[:,:,14],admmst.s[:,:,14],admmst.taus[:,:,14],a[14])
freq15 = genshared3D(admmst.x[:,:,15],admmst.p[:,:,15],admmst.taup[:,:,15],admmst.wlt[:,:,15],psfst.mypsf[:,:,15],admmst.fty[:,:,15],admmst.s[:,:,15],admmst.taus[:,:,15],a[15])

freq1_tt = genshared4D(admmst.t[:,:,1,:],admmst.taut[:,:,1,:],nspat,a[1])
freq2_tt = genshared4D(admmst.t[:,:,2,:],admmst.taut[:,:,2,:],nspat,a[2])
freq3_tt = genshared4D(admmst.t[:,:,3,:],admmst.taut[:,:,3,:],nspat,a[3])
freq4_tt = genshared4D(admmst.t[:,:,4,:],admmst.taut[:,:,4,:],nspat,a[4])
freq5_tt = genshared4D(admmst.t[:,:,5,:],admmst.taut[:,:,5,:],nspat,a[5])
freq6_tt = genshared4D(admmst.t[:,:,6,:],admmst.taut[:,:,6,:],nspat,a[6])
freq7_tt = genshared4D(admmst.t[:,:,7,:],admmst.taut[:,:,7,:],nspat,a[7])
freq8_tt = genshared4D(admmst.t[:,:,8,:],admmst.taut[:,:,8,:],nspat,a[8])
freq9_tt = genshared4D(admmst.t[:,:,9,:],admmst.taut[:,:,9,:],nspat,a[9])
freq10_tt = genshared4D(admmst.t[:,:,10,:],admmst.taut[:,:,10,:],nspat,a[10])
freq11_tt = genshared4D(admmst.t[:,:,11,:],admmst.taut[:,:,11,:],nspat,a[11])
freq12_tt = genshared4D(admmst.t[:,:,12,:],admmst.taut[:,:,12,:],nspat,a[12])
freq13_tt = genshared4D(admmst.t[:,:,13,:],admmst.taut[:,:,13,:],nspat,a[13])
freq14_tt = genshared4D(admmst.t[:,:,14,:],admmst.taut[:,:,14,:],nspat,a[14])
freq15_tt = genshared4D(admmst.t[:,:,15,:],admmst.taut[:,:,15,:],nspat,a[15])


#########################################################################
# genshared3D(x,p,taup,wlt,psf,fty,s,taus)



function genshared3D(x::Array{Float64,3},p::Array{Float64,3},taup::Array{Float64,3},wlt::Array{Float64,3},psf::Array{Float64,3},
                     fty::Array{Float64,3},s::Array{Float64,3},taus::Array{Float64,3})


    Ndim = 8
    nxy = size(x)[1]
    nfreq = size(x)[3]

    workers = nworkers()
    a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
    a = a[1:nfreq]


    tabname = "freq"
    toeval = ""

    listarr = Dict{Any,Any}
    listarr = {0 => 1}
    for z in 1:nfreq
          listarr[z] = SharedArray(Float64,nxy,nxy,Ndim,pids=[1,a[z]]);
    end



    for z in 1:nfreq
        toeval = string(toeval,string(tabname, "$z", "=", "similar(listarr[$z]);"))
        println(toeval)
    end
    println(" 1")
    println(listarr[1])
    println("2 ")
    println(listarr[2])
    println(" 3")
    println(listarr[3])
    println(" 4")
    println(listarr[4])
    println(parse(toeval))
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

        chaine = string("for nd in 1:Ndim;",string(tabname, "$z", "[:,:,nd]","=","$result[nd];")," end;")

        toeval = string(toeval,chaine)
    end
    eval(parse(toeval))

end


#########################################################################
Ndim = 8
nxy = 2
nfreq = 4

workers = nworkers()
a = int(repmat(linspace(2,workers+1,workers),int(ceil(nfreq/workers)),1))
a = a[1:nfreq]

x = randn(nxy,nxy,nfreq);
p =randn(nxy,nxy,nfreq);
taup = randn(nxy,nxy,nfreq);
wlt = randn(nxy,nxy,nfreq);
psf =randn(nxy,nxy,nfreq);
fty = randn(nxy,nxy,nfreq);
s =randn(nxy,nxy,nfreq);
taus =randn(nxy,nxy,nfreq);


tabname = "freq"
toeval = ""

listarr = {0 => 1}
for z in 1:nfreq
      listarr[z] = SharedArray(Float64,nxy,nxy,Ndim,pids=[1,a[z]]);
end

for n in 1:nfreq
    proc = n
  toeval = string(toeval,string(tabname, "$n", "=", "similar(listarr[$n]);"))
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

    chaine = string("for nd in 1:Ndim;",string(tabname, "$z", "[:,:,nd]","=","$result[nd];")," end;")

    toeval = string(toeval,chaine)
end
eval(parse(toeval))
#########################################################################
