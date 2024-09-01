using Combinatorics
using Distributions
using SpecialFunctions








function padArray(array::Vector{},len::Int)::Vector{Int64}
    return(append!(array,zeros(len-size(array)[1])))
end
function mirrorList(oldList, length::Int)
    newList = [oldList; reverse(oldList)]
    if (length%2==1)
        newList = deleteat!(newList,trunc(Int,length÷2)+1)
    end
    return newList
end
function getProbabilities(ham,time)
    T = exp(-im*time*ham)
    probabilities = (real(T).^2)+(imag(T).^2)
    return probabilities
end
function randBoundedFill(bounds)
    arr = []
    for i in 1:length(bounds[1])
        append!(arr,rand(Uniform(bounds[1][i],bounds[2][i])))
    end
    return arr
end









function makeCavBasis(nexc::Int,ncav::Int)::Vector{Vector{Int8}}
    if nexc==0
        return [zeros(ncav)]
    end
    basis=integer_partitions(nexc) # There are n total particles. Divide them up into possble permutations
    basis=filter(x->size(x)[1]<=ncav,basis) # remove paritions with more necessary locations than cavities we have
    basis=map(x -> append!(x,zeros(ncav-size(x)[1])),basis) # add zeros to get number of cavities
    total_basis=[]
    for i in basis
        append!(total_basis,multiset_permutations(i,ncav)) # add all rearrangements of possible configurations
    end
    return total_basis
end
function makeTotalBasis(nexc::Int,ncav::Int,nemi::Int)
    total_basis=[] # important to note: emitters can only contain either 0 or 1 excitations.
    for i in 0:min(nexc,nemi) # until running out of emitters or excitations
        emibasis=multiset_permutations(padArray(ones(i),nemi),nemi) # get ways emitters can be excited
        cavbasis=makeCavBasis(nexc-i,ncav) # and the cavity basis for the remaining excitations
        for emi_i in emibasis
            for cav_j in cavbasis
                append!(total_basis,[[cav_j;emi_i]]) # add all combinations of emitter and cavity congifurations to basis
            end
        end
    end
    return(total_basis)
end
function makeHam2(cav_J,emi_g,n_exc; emi_idx=[])
    n_cav=size(cav_J)[1]+1
    n_emi=size(emi_g)[1]
    if emi_idx == []
        emi_index = collect(1:n_emi) # if optional element emi_idx is left out assume 1 emitter in all cavites
    else
        emi_index = emi_idx
    end
    basis=makeTotalBasis(n_exc,n_cav,n_emi) # start with a basis
    n_basis=size(basis)[1]
    ham=zeros(n_basis,n_basis) # initialize a Ham filled with 0s
    for c_basis in eachindex(basis) # for each state in the basis
        indexs=findall(!iszero, basis[c_basis][1:n_cav]) 
        for c_index in  indexs # for every cavity that has an excitation
            if c_index>1 
                tempbasis=copy(basis[c_basis]) # make a new state with the excitation moved one to the left
                tempbasis[c_index]-=1
                tempbasis[c_index-1] +=1
                cpoint=findfirst(==(tempbasis),basis) # and find the index of that state
                ham[c_basis,cpoint]+=cav_J[c_index-1]*sqrt(basis[c_basis][c_index])*sqrt(tempbasis[c_index-1])
                # connect these states in the Ham with the appropriate J for the hop * the relevant coefficents
            end
            if c_index<n_cav # do same connection for a hop to the right
                tempbasis=copy(basis[c_basis])
                tempbasis[c_index]-=1
                tempbasis[c_index+1] +=1
                cpoint=findfirst(==(tempbasis),basis)
                ham[c_basis,cpoint]+=cav_J[c_index]*sqrt(basis[c_basis][c_index])*sqrt(tempbasis[c_index+1])
            end
        end
    end
    for c_basis in eachindex(basis) # add g connection for excitation entering  emitter 
        for i in 1:n_emi
            if basis[c_basis][n_cav+i]==0 #if emitter is empty
                if basis[c_basis][emi_index[i]]>0 # and cavity has excitation
                    tempbasis=copy(basis[c_basis])
                    tempbasis[emi_index[i]]-=1
                    tempbasis[n_cav+i]+=1
                    cpoint=findfirst(==(tempbasis),basis) # find state after excitation
                    ham[c_basis,cpoint]+=sqrt(basis[c_basis][emi_index[i]])*emi_g[i] # connect with g * coefficent
                end
            end
        end
        for i in 1:n_emi # same process but for exitation leaving emitter
            if basis[c_basis][n_cav+i]==1
                tempbasis=copy(basis[c_basis])
                tempbasis[emi_index[i]]+=1
                tempbasis[n_cav+i]-=1
                cpoint=findfirst(==(tempbasis),basis)
                ham[c_basis,cpoint]+=sqrt(tempbasis[emi_index[i]])*emi_g[i]
            end
        end
    end
    return(ham)
end

function makeSpinlessBasis(n_sites; n_exc=1) # makes permutations of possible electron configurations
    basis = [fill(0.,2*n_sites)]
    for n in 1:(2*n_sites*n_exc)
        partitions = integer_partitions(n)
        legalPartitions = Vector{Float64}[]
        for p in partitions
            if !(size(p)[1]>2*n_sites || any(x->x>n_exc,p))
                append!(legalPartitions, multiset_permutations(padArray(p,2*n_sites),2*n_sites))
            end
        end
        append!(basis, legalPartitions)
    end
    return basis
end
function getTotalElectrons(state)
    n_sites = Int(size(state)[1]/2)
    total = 0
    for i in 1:2*n_sites
        total += state[i]
    end
    return total
end

function filterBasis(basis, nup, ndown)
    newBasis = Array{Float64}[]
    n_sites = trunc(Int64,size(basis[1])[1]/2)
    for state in basis
        if getTotalElectrons(state[1:n_sites])!=nup || getTotalElectrons(state[n_sites+1:2*n_sites])!=ndown 
            continue
        end
        append!(newBasis, [state])
    end
    return newBasis
end
function makeAnyHam(name,n_sites,n_up,n_down,t,v,U;display=false, returnBasis=false, basis=[]) 
    if name in ["Hubbard", "Hub", "hub", "H", "h"]
        haveTwoRows,linkRowsVertical,linkBottomRow = false,false,false
        horizontalHops = vcat(1:n_sites-1,n_sites:2*n_sites-1)
    elseif name in ["Periodic Anderson Model", "PAM", "Pam", "pam", "P", "p"] # Pam upper U should be 0
        haveTwoRows,linkRowsVertical,linkBottomRow = true,true,false
        horizontalHops = vcat(1:n_sites-1,2*n_sites:3*n_sites-1)
    elseif name in ["Ladder", "Lad", "lad", "L", "l"]
        haveTwoRows,linkRowsVertical,linkBottomRow = true,true,true
        horizontalHops = vcat(1:n_sites-1,n_sites:2*n_sites-1,2*n_sites:3*n_sites-1,3*n_sites:4*n_sites-1)
    elseif name in ["Jaynes-Cummings-Hubbard", "JCH", "JCHH", "J", "j"]
        return makeHam2(t,v,n_up)
    else
        throw(DomainError("Hamiltonian type not found"))
    end
    if basis == []
        if haveTwoRows
            basis = makeSpinlessBasis(n_sites*2)
        else
            basis = makeSpinlessBasis(n_sites)
        end
        basis = filterBasis(basis,n_up,n_down)
    end
    if display
        showBasis(basis,name=name)
    end
    n_states = size(basis)[1]
    ham = fill(0.0, (n_states,n_states))

    for s in 1:n_states
        state = basis[s]
        
        # Diagonal elements from U interaction
        diag = 0.0
        if haveTwoRows
            for i in vcat(1:2*n_sites)
                if state[i]==1.0 && state[i+(2*n_sites)]==1.0
                     diag += U[(i-1)÷(n_sites)+1]
                end
            end
        else
            for i in 1:n_sites
                if state[i]==1.0 && state[i+n_sites]==1.0
                    diag += U[1]
                end
            end
        end

        ham[s,s] = diag

        # Do allowed horizontal hops
        for i in horizontalHops
            if i % n_sites != 0 && state[i] == 1 && state[i+1] == 0 # electron hop right
                tempstate = copy(state)
                tempstate[i+1] = 1.0
                tempstate[i] = 0.0
                idx = findfirst(x->x==tempstate,basis)
                if idx !=nothing
                    ham[s,idx] = t[i%n_sites]
                    ham[idx,s] = t[i%n_sites]
                end
            end
        end

        # link top and bottoms
        if linkRowsVertical
            for i in vcat(1:n_sites,2*n_sites+1:3*n_sites)
                if ((i-1) % (2*n_sites) < n_sites) && state[i] == 1 && state[i+n_sites] == 0 # electron hop down
                    tempstate = copy(state)
                    tempstate[i+n_sites] = 1.0
                    tempstate[i] = 0.0
                    count = 0
                    for j in i+1:i+n_sites-1 # look between the sites
                        if state[j] == 1.0 # count the up electrons we 'hop over'
                            count += 1
                        end
                    end
                    idx = findfirst(x->x==tempstate,basis)
                    if idx !=nothing
                        ham[s,idx] = v[(i-1) % n_sites + 1] * ((-1)^count) # here are the 'second quantization sign flips'
                        ham[idx,s] = v[(i-1) % n_sites + 1] * ((-1)^count)
                    end
                end
            end
        end
        
    end
    if returnBasis
        return ham, basis
    else
        return ham
    end
end








function randBoundedFill(bounds)
    arr = []
    for i in 1:length(bounds[1])
        append!(arr,rand(Uniform(bounds[1][i],bounds[2][i])))
    end
    return arr
end
function dynamicOutput(txtout,outstring)
    if isnothing(txtout)
        println(outstring)
    else
        file = open(txtout, "a")
        write(file, (outstring))
        close(file)
        println(outstring)
    end
end









function dualAnnealing(numVars::Int64, errFunc::Function, maxIterations::Int64; initialTemp=5230.0,
        tempRestartRatio=2.e-5, visitingBias=2.62, acceptanceBias=-5.0, maxIterSinceImprovement=1000,
        seed=nothing, bounds=(nothing,nothing), start=nothing, evaluateSolution::Function=identity, errorBenchmarks=[],
        localSearches=0, printIterations=0, txtout=nothing, seperateRuns=false, graphConvergence=false, args...)
    # invalid input detection
    bounds==(nothing,nothing) && isnothing(start) && throw(DomainError("No start or bounds have been provided."))
    !(length(bounds[1])==length(bounds[2])) && throw(DomainError("Lower and upper bounds are unequal sizes"))
    !(all(bounds[1].<bounds[2])) && throw(DomainError("Lower bounds are not all lower than upper bounds"))
    !(initialTemp > 0) && throw(DomainError("Temperature is not positive"))
    !(0 < tempRestartRatio < 1) && throw(DomainError("tempRestartRatio is not between 0 and 1"))
    !(1 < visitingBias <= 3) && throw(DomainError("visitingBias is not in the set (1, 3]"))
    !(-10000 < acceptanceBias <= -5) && throw(DomainError("acceptanceBias is not in the set (-10000, -5]"))
    !(isnothing(start)) && !(length(start)==numVars)
    
    boundsRange = bounds[2] .- bounds[1]
    evaluatorProvided = !(evaluateSolution==identity)
    if graphConvergence
        recordedIterations = round.(collect(0:0.005:1).*maxIterations)
        convergenceCurrentData = []
        convergenceBestData = []
        convergenceTemperatureData = []
    end
    if !(isnothing(txtout)) && seperateRuns
        file = open(txtout, "a")
        seperateRuns && write(file,"STARTING NEW RUN\n")
        close(file)
    end
    !(isnothing(seed)) && Random.seed!(seed)
    if isnothing(start)
        currentLoc = randBoundedFill(bounds)
    else
        currentLoc = start
    end
    currentErr = errFunc(currentLoc; args...)
    bestLoc,bestErr = currentLoc,currentErr
    tempRestartThreshold = initialTemp * tempRestartRatio
    iteration = 0
    tempIteration = 0
    iterSinceImprovement = 0
    runLoop = true
    t1 = exp((visitingBias - 1) * log(2.0)) - 1.0
    functionCalls = 1
    
    # Main Optimization Loop
    while(runLoop)
        iteration >= maxIterations && (runLoop = false)
        # Temperature update
        t2 = exp((visitingBias - 1) * log(tempIteration + 2.0)) - 1.0
        temp = initialTemp * t1 / t2
        # Printing updates
        if !(printIterations==0) && iteration%printIterations==0
            outstring = "Iteration: "*string(iteration)*". Best Error: "*string(round(bestErr,digits=8))*
                ". Current Error: "*string(round(currentErr,digits=8))*". Current Location: "*
                string(round.(currentLoc,digits=8))*"\n"
            dynamicOutput(txtout,outstring)
            GC.gc()
        end
        # Convergence data recording
        if graphConvergence && (iteration in recordedIterations)
            append!(convergenceCurrentData, currentErr)
            append!(convergenceBestData, bestErr)
            append!(convergenceTemperatureData,temp)
        end
        # Reannealing process. In the python Dual Annealing tmeperature is not reset. I believe this is a mistake
        if temp < tempRestartThreshold
#             !(printIterations==0) && println("Temperature and location are reset") # Make depend on printInfo variable
            tempIteration = 0
            if isnothing(start) 
                currentLoc = randBoundedFill(bounds)
                currentErr = errFunc(currentLoc; args...)
            end
        end
        subTemp = temp / (iteration + 1)

        # Create visiting distribution based on temperature. Cauchy-Lorentz distribution
        # [2] Tsallis C, Stariolo DA. Generalized Simulated Annealing.
        # Physica A, 233, 395-406 (1996).
        factor1 = exp(log(temp) / (visitingBias - 1.0))
        factor2 = exp((4.0 - visitingBias) * log(visitingBias - 1.0))
        factor3 = exp((2.0 - visitingBias) * log(2.0) / (visitingBias - 1.0))
        factor4p = sqrt(pi) * factor2 / (factor3 * (3.0 - visitingBias))
        factor4 = factor4p * factor1
        factor5 = 1.0 / (visitingBias - 1.0) - 0.5
        factor6 = pi * (1.0 - factor5) / sin(pi * (1.0 - factor5)) / exp(logabsgamma(2.0 - factor5)[1])
        
        # Strategy chain
        iterSinceImprovement += 1
        chainMinLoc, chainMinErr = currentLoc, currentErr
        for j in 1:numVars*2
            # Create a random visit based on our distribution
            x, y = rand(Normal(),numVars),rand(Normal(),numVars)
            x *= exp(-(visitingBias - 1.0) * log(factor6 / factor4) / (3.0 - visitingBias))
            den = exp.((visitingBias - 1.0) * log.(abs.(y)) / (3.0 - visitingBias))
            visits = x ./ den
            visits = map(x-> abs(x) > 1e8 ? sign(x)*1e8*rand() : x, visits)
            # Change all variables or an individual variable based on case
            if j <= numVars
                visitLoc = currentLoc .+ visits
            else
                visitLoc = deepcopy(currentLoc)
                visitLoc[(j%numVars) + 1] = visits[(j%numVars) + 1]
            end
            # Keep visits in bounds. Going over one bound wraps to the other
            if !(bounds==(nothing,nothing))
                a = visitLoc .- bounds[1]
                b = (a .% boundsRange) .+ boundsRange
                visitLoc = (b .% boundsRange) .+ bounds[1]
            end
            visitErr = errFunc(visitLoc; args...)
            functionCalls += 1
            # Acceptance if visit is better
            if visitErr < currentErr 
                currentLoc, currentErr = visitLoc, visitErr
                if visitErr < bestErr
                   for benchmark in errorBenchmarks
                        if !isnothing(errorBenchmarks) && benchmark < bestErr && benchmark >= visitErr
                            outstring = "Iteration: "*string(iteration)*". Function Calls: "*string(functionCalls)*
                                ". Reached Error Benchmark: "*string(benchmark)*". Current Error: "*
                                string(round(currentErr,digits=8))*". Current Location: "*
                                string(round.(currentLoc,digits=8))*"\n"
                            dynamicOutput(txtout,outstring)
                        end
                    end
                    # Accept Change
                    bestLoc, bestErr = visitLoc, visitErr
                    iterSinceImprovement = 0
                    doLocalSearch = true
                    # If new solution meets the provided evaluation we can stop the loop 
                    evaluatorProvided && (runLoop = !evaluateSolution(bestLoc; args...))

                end
            # Random acceptance or rejection if visit is worse
            else
                decisionRatio = 1.0 - ((1.0 - acceptanceBias) * (visitErr - currentErr) / subTemp)
                if decisionRatio < 0
                    decisionRatio = 0
                else
                    decisionRatio = exp(log(decisionRatio) / (1.0 - acceptanceBias))
                end
                randDecision = rand()
                if randDecision <= decisionRatio
                    currentLoc, currentErr = visitLoc, visitErr
                end
            end
            # If no new minimum has been found in many iterations, find strategy chain minimum for local search
            if iterSinceImprovement >= maxIterSinceImprovement
                visitErr < chainMinErr && (chainMinLoc, chainMinErr = visitLoc, visitErr)
            end
        end
        
        # Local search
        if !(localSearches==0)
            doLocalSearch = false
            iterSinceImprovement >= maxIterSinceImprovement && (doLocalSearch = true)
            # Chance of doing local search anyways
            localSearchChance = exp(100 * numVars * (bestErr - currentErr) / subTemp)
            rand() < localSearchChance && (doLocalSearch=true)
            if doLocalSearch
                #Add in local search
                # "Local search" is actually a gradient descent (actually L-BFGS-B gradient approximation)
                println("doing a local search")
            end
        end
        
        tempIteration += 1
        iteration += 1
    end
    if !(printIterations==0)
        outstring = "Iteration: "*string(iteration)*". Best Error: "*string(round(bestErr,digits=8))*
            ". Current Error: "*string(round(currentErr,digits=8))*". Best Location: "*
            string(round.(bestLoc,digits=8))*"\n"
        dynamicOutput(txtout,outstring)
    end
    if graphConvergence
        recordedIterations = recordedIterations[1:length(convergenceCurrentData)]
        display(plot(recordedIterations,[convergenceCurrentData,convergenceBestData,convergenceTemperatureData],
                yaxis=:log, label=["Current Error" "Best Error" "Temperature"]))
    end
    
    return bestLoc, bestErr
end












function OneInputSingleJCHHTransferError(x; cavityLength=0, startIdx=1, goalIdx=1, n_exc=1)
    cavityLength==0 && throw(DomainError("The system cavity length must be specified"))
    Jnumber = trunc(Int,(cavityLength)*0.5)
    gnumber = trunc(Int,(cavityLength+1)*0.5)
    myJ = mirrorList(x[1:Jnumber],cavityLength-1)
    myg = mirrorList(x[Jnumber+1:Jnumber+gnumber],cavityLength)
    time = x[Jnumber+gnumber+1]
    fullHam = makeHam2(myJ,myg,n_exc)
    error = 1.0 - getProbabilities(fullHam,time)[startIdx,goalIdx]
    minimum(x) < (sum(x)/length(x))/10 && (error = error+1.0) # Make very low values 'invalid'
    maximum(x) > (sum(x)/length(x))*10 && (error = error+1.0) # Make very high values 'invalid'
    return error
end
function myErrorCallback(x; cavityLength=0, startIdx=1, goalIdx=1, n_exc=1)
    err = OneInputSingleJCHHTransferError(x, cavityLength=cavityLength, startIdx=startIdx, goalIdx=goalIdx, n_exc=n_exc)
    err > 0.01 && return false
    return true
end
function myErrorCallbackMultipleTransfers(x; showResult=false, transfers=2)
    err = JCHH8ErrorMultipleTransfers(x, showResult=showResult, transfers=transfers)
    err > 0.01 && return false
    return true
end
function myErrorCallbackMultipleTransfersReduced(x; showResult=false, transfers=2)
    err = JCHH6ErrorMultipleTransfers(x, showResult=showResult, transfers=transfers)
    err > 0.02 && return false
    return true
end
function OneInputSinglePAMTransferError(x; cavityLength=0, startIdx=1, goalIdx=1, electrons=(1,1))
    cavityLength==0 && throw(DomainError("The system cavity length must be specified"))
    Jnumber = trunc(Int,(cavityLength)*0.5)
    gnumber = trunc(Int,(cavityLength+1)*0.5)
    myt = mirrorList(x[1:Jnumber],cavityLength-1)
    myv = mirrorList(x[Jnumber+1:Jnumber+gnumber],cavityLength)
    myU = x[Jnumber+gnumber+1:Jnumber+gnumber+2]
    time = x[Jnumber+gnumber+3]
    fullHam = makeAnyHam("PAM",cavityLength,electrons[1],electrons[2],myt,myv,myU)
    error = 1.0 - getProbabilities(fullHam,time)[startIdx,goalIdx]
    minimum(x) < (sum(x)/length(x))/10 && (error = error+1.0) # Make very low values 'invalid'
    maximum(x) > (sum(x)/length(x))*10 && (error = error+1.0) # Make very high values 'invalid'
    return error
end
function myErrorCallbackPAM(x; cavityLength=0, startIdx=1, goalIdx=1, electrons=(1,1))
    err = OneInputSinglePAMTransferError(x, cavityLength=cavityLength, startIdx=startIdx, goalIdx=goalIdx, electrons=electrons)
    return err < 0.01
end


function JCHH8ErrorMultipleTransfers(x; showResult=false, transfers=1)
    TransferList = [[29,36],[37,100],[1,28],[38,99]] #Fill this with 8 pairs
    myt = mirrorList(x[1:4],7)
    myv = mirrorList(x[5:8],8)
    myU = []
    times = x[9:8+transfers]
    fullHam = makeHam2(myt,myv,2; emi_idx=[])
    error = 0.0
    for i in 1:transfers
        error +=1.0 - getProbabilities(fullHam,times[1])[TransferList[i][1],TransferList[i][2]]
        showResult && println("transfer from "*string(TransferList[i][1])*" to "*string(TransferList[i][2])*
            " has fidelity: "*string(1.0-error))
    end
    minimum(x) < (sum(x)/length(x))/10 && (error = error+1.0) # Make very low values 'invalid'
    maximum(x) > (sum(x)/length(x))*10 && (error = error+1.0) # Make very high values 'invalid'
    return error
end
function JCHH6ErrorMultipleTransfers(x; showResult=false, transfers=2)
    TransferList = [[16,21],[22,57],[1,15],[58,72]]
    myt = mirrorList(x[1:3],5)
    myv = mirrorList(x[4:6],6)
    myU = []
    times = x[7:6+transfers]
    fullHam = makeHam2(myt,myv,2; emi_idx=[])
    error = 0.0
    for i in 1:transfers
        error +=1.0 - getProbabilities(fullHam,times[1])[TransferList[i][1],TransferList[i][2]]
        showResult && println("transfer from "*string(TransferList[i][1])*" to "*string(TransferList[i][2])*
            " has fidelity: "*string(1.0-error))
    end
    error = error / transfers
    minimum(x) < (sum(x)/length(x))/10 && (error = error+1.0) # Make very low values 'invalid'
    maximum(x) > (sum(x)/length(x))*10 && (error = error+1.0) # Make very high values 'invalid'
    return error
    # 16  |2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>  
    # 21  |0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0> 
    
    # 22  |1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0>
    # 57  |0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1> 

    #  1  |1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>
    # 15  |0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0>

    # 58  |0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0> 
    # 72  |0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1>  
end
function JCHH4ErrorMultipleTransfers(x; showResult=false, transfers=2)
    TransferList = [[7,10],[11,26],[1,6],[27,32]]
    myt = mirrorList(x[1:2],3)
    myv = mirrorList(x[3:4],4)
    myU = []
    times = x[5:4+transfers]
    fullHam = makeHam2(myt,myv,2; emi_idx=[])
    error = 0.0
    for i in 1:transfers
        error +=1.0 - getProbabilities(fullHam,times[1])[TransferList[i][1],TransferList[i][2]]
        showResult && println("transfer from "*string(TransferList[i][1])*" to "*string(TransferList[i][2])*
            " has fidelity: "*string(1.0-error))
    end
    error = error / transfers
    minimum(x) < (sum(x)/length(x))/10 && (error = error+1.0) # Make very low values 'invalid'
    maximum(x) > (sum(x)/length(x))*10 && (error = error+1.0) # Make very high values 'invalid'
    return error
    # 7  |2, 0, 0, 0, 0, 0, 0, 0> 
    # 10  |0, 0, 0, 2, 0, 0, 0, 0>
    
    # 11  |1, 0, 0, 0, 1, 0, 0, 0> 
    # 26  |0, 0, 0, 1, 0, 0, 0, 1> 

    #  1  |1, 1, 0, 0, 0, 0, 0, 0> 
    #  6  |0, 0, 1, 1, 0, 0, 0, 0> 

    # 27  |0, 0, 0, 0, 1, 1, 0, 0> 
    # 32  |0, 0, 0, 0, 0, 0, 1, 1>
end
function JCHH5ErrorMultipleTransfers(x; showResult=false, transfers=1)
    TransferList = [[1,10],[2,9],[3,7],[5,8],[11,15],[12,14]]
    myt = mirrorList(x[1:2],4)
    myv = mirrorList(x[3:5],5)
    myU = []
    times = x[6:5+transfers]
    fullHam = makeHam2(myt,myv,2; emi_idx=[])
    error = 0.0
    for i in 1:transfers
        error +=1.0 - getProbabilities(fullHam,times[1])[TransferList[i][1],TransferList[i][2]]
        showResult && println("transfer from "*string(TransferList[i][1])*" to "*string(TransferList[i][2])*
            " has fidelity: "*string(1.0-error))
    end
    error = error / transfers
    minimum(x) < (sum(x)/length(x))/10 && (error = error+1.0) # Make very low values 'invalid'
    maximum(x) > (sum(x)/length(x))*10 && (error = error+1.0) # Make very high values 'invalid'
    return error
    # 1  |1, 1, 0, 0, 0, 0, 0, 0, 0, 0> 
    #10  |0, 0, 0, 1, 1, 0, 0, 0, 0, 0> 

    # 2  |1, 0, 1, 0, 0, 0, 0, 0, 0, 0> 
    # 9  |0, 0, 1, 0, 1, 0, 0, 0, 0, 0> 

    # 3  |1, 0, 0, 1, 0, 0, 0, 0, 0, 0>
    # 7  |0, 1, 0, 0, 1, 0, 0, 0, 0, 0> 

    # 5  |0, 1, 1, 0, 0, 0, 0, 0, 0, 0> 
    # 8  |0, 0, 1, 1, 0, 0, 0, 0, 0, 0> 

    # 11  |2, 0, 0, 0, 0, 0, 0, 0, 0, 0> 
    # 15  |0, 0, 0, 0, 2, 0, 0, 0, 0, 0> 

    # 12  |0, 2, 0, 0, 0, 0, 0, 0, 0, 0> 
    # 14  |0, 0, 0, 2, 0, 0, 0, 0, 0, 0>
end
function JCHH4ErrorMultipleTransfersRandom(x; showResult=false, transfers=1)
    TransferList = [[1, 6], [2, 5], [7, 10], [8, 9], [11, 26], [12, 25], [13, 24], [14, 23], [15, 22], [16, 21], 
        [17, 20], [18, 19], [27, 32], [28, 31]] #Fill this with pairs
    myRandTransfers = [[-1,-1] for i in 1:transfers]
    sample!(TransferList, myRandTransfers; replace=false)
    myt = mirrorList(x[1:2],3)
    myv = mirrorList(x[3:4],4)
    myU = []
    times = x[5:4+transfers]
    fullHam = makeHam2(myt,myv,2; emi_idx=[])
    error = 0.0
    for i in 1:transfers
        error +=1.0 - getProbabilities(fullHam,times[1])[myRandTransfers[i][1],myRandTransfers[i][2]]
        showResult && println("transfer from "*string(myRandTransfers[i][1])*" to "*string(myRandTransfers[i][2])*
            " has fidelity: "*string(1.0-error))
    end
    minimum(x) < (sum(x)/length(x))/10 && (error = error+1.0) # Make very low values 'invalid'
    maximum(x) > (sum(x)/length(x))*10 && (error = error+1.0) # Make very high values 'invalid'
    return error
end

function JCHH4ErrorMultipleTransfersRandomSimultaneous(x; showResult=false, transfers=1)
    TransferList = [[1, 6], [2, 5], [7, 10], [8, 9], [11, 26], [12, 25], [13, 24], [14, 23], [15, 22], [16, 21], 
        [17, 20], [18, 19], [27, 32], [28, 31]] #Fill this with pairs
    myRandTransfers = [[-1,-1] for i in 1:transfers]
    sample!(TransferList, myRandTransfers; replace=false)
    myt = mirrorList(x[1:2],3)
    myv = mirrorList(x[3:4],4)
    myU = []
    time = x[5]
    fullHam = makeHam2(myt,myv,2; emi_idx=[])
    error = 0.0
    for i in 1:transfers
        error +=1.0 - getProbabilities(fullHam,time)[myRandTransfers[i][1],myRandTransfers[i][2]]
        showResult && println("transfer from "*string(myRandTransfers[i][1])*" to "*string(myRandTransfers[i][2])*
            " has fidelity: "*string(1.0-error))
    end
    minimum(x) < (sum(x)/length(x))/10 && (error = error+1.0) # Make very low values 'invalid'
    maximum(x) > (sum(x)/length(x))*10 && (error = error+1.0) # Make very high values 'invalid'
    return error
end





# for numcav in [3,4,5]
#     randomString = string(trunc(Int,rand()*9000+1000))
#     mytextout="JCHH"*string(numcav)*"SST_ConvergenceV3"*randomString*".txt"
#     for i in 1:8
#         bestX, bestXerr = dualAnnealing(numcav+1,OneInputSingleJCHHTransferError,200000,bounds=(fill(0.,numcav+1), fill(300.,numcav+1)),
#             printIterations=1000,graphConvergence=false,seperateRuns=true, txtout=mytextout, errorBenchmarks=[0.05,0.02,0.01],
#             evaluateSolution=myErrorCallback, cavityLength=numcav, startIdx=Int(numcav*(numcav-1)/2 + 1), goalIdx=Int(numcav*(numcav+1)/2), n_exc=2)
#         println("Run number: ",i," The best solution was: \n", round.(bestX,digits=8), "\nWith error: ",bestXerr)
#     end
# end

# randomString = string(trunc(Int,rand()*9000+1000))
# mytextout="JCHH4_2ST_HighFiData"*randomString*".txt"
# for i in 1:20
#     bestX, bestXerr = dualAnnealing(6,JCHH4ErrorMultipleTransfers,50000,bounds=(fill(0.,6), fill(300.,6)),
#         printIterations=1000,graphConvergence=false,seperateRuns=true, txtout=mytextout, errorBenchmarks=[0.05,0.02,0.01,0.001],
#         transfers=2)
#     println("Run number: ",i," The best solution was: \n", round.(bestX,digits=8), "\nWith error: ",bestXerr)
# end



numCav = 5
transferList = [[1, numCav*(2*numCav - 1)], [numCav, 2*numCav*(numCav-1) + 1]]

for i in 2:2
    transfer = i
    mytextout="PAM"*string(numCav)*"_"*string(transfer)*".txt"
    println(mytextout)

    bestX, bestXerr = dualAnnealing(numCav+3,OneInputSinglePAMTransferError,50000,bounds=(fill(0.,numCav+3), fill(300.,numCav+3)),
        printIterations=1000,graphConvergence=false,seperateRuns=true, txtout=mytextout, errorBenchmarks=[0.05,0.02,0.01],
        evaluateSolution=myErrorCallbackPAM, cavityLength=numCav, startIdx=transferList[transfer][1], goalIdx=transferList[transfer][2], electrons=(1,1))
    println("Run number: ",i," The best solution was: \n", round.(bestX,digits=8), "\nWith error: ",bestXerr)
end



# randomString = string(trunc(Int,rand()*9000+1000))
# mytr = 4
# mytextout="JCHH6SST_TS4_ConvergenceV2"*randomString*".txt"
# println(mytextout)
# for i in 1:15
#     bestX, bestXerr = dualAnnealing(6+mytr,JCHH6ErrorMultipleTransfers,50000,bounds=(fill(0.,6+mytr), fill(300.,6+mytr)),# 
#         printIterations=1000,graphConvergence=false,seperateRuns=true, txtout=mytextout, errorBenchmarks=[0.1,0.05,0.02],
#         evaluateSolution=myErrorCallbackMultipleTransfersReduced, transfers=mytr)
#     println("Run number: ",i," The best solution was: \n", round.(bestX,digits=8), "\nWith error: ",bestXerr)
# end


# for mytr in [1,3,4,5,6]
#     randomString = string(trunc(Int,rand()*9000+1000))
#     mytextout="JCHH5_"*string(mytr)*"transfers"*randomString*".txt"
#     println(mytextout)
#     for i in 1:5
#         bestX, bestXerr = dualAnnealing(5+mytr,JCHH5ErrorMultipleTransfers,50000,bounds=(fill(0.,5+mytr), fill(300.,5+mytr)),
#             printIterations=1000,graphConvergence=false,seperateRuns=true, txtout=mytextout, errorBenchmarks=[0.05,0.02,0.01,0.001],
#              transfers=mytr)
#         println("Run number: ",i," The best solution was: \n", round.(bestX,digits=8), "\nWith error: ",bestXerr)
#     end
# end

# for mytr in [1,2,3,4]
#     randomString = string(trunc(Int,rand()*9000+1000))
#     mytextout="JCHH8_"*string(mytr)*"transfers"*randomString*".txt"
#     println(mytextout)
#     for i in 1:2
#         bestX, bestXerr = dualAnnealing(8+mytr,JCHH8ErrorMultipleTransfers,50000,bounds=(fill(0.,8+mytr), fill(300.,8+mytr)),
#             printIterations=1000,graphConvergence=false,seperateRuns=true, txtout=mytextout, errorBenchmarks=[0.05,0.02,0.01,0.001],
#              transfers=mytr)
#         println("Run number: ",i," The best solution was: \n", round.(bestX,digits=8), "\nWith error: ",bestXerr)
#     end
# end

# for mytr in [2,3,4,5,6]
#     randomString = string(trunc(Int,rand()*9000+1000))
#     mytextout="JCHH4_"*string(mytr)*"randomTransfers"*randomString*".txt"
#     println(mytextout)
#     for i in 1:10
#         bestX, bestXerr = dualAnnealing(4+mytr,JCHH4ErrorMultipleTransfersRandom,50000,bounds=(fill(0.,4+mytr), fill(300.,4+mytr)),
#             printIterations=1000,graphConvergence=false,seperateRuns=true, txtout=mytextout, errorBenchmarks=[0.05,0.02,0.01,0.001],
#              transfers=mytr)
#         println("Run number: ",i," The best solution was: \n", round.(bestX,digits=8), "\nWith error: ",bestXerr)
#     end
# end