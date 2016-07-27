module LCIOUtils
export getpositions, mappositions, mapparticles, maptruthparticles, numfaces, hcal_outer, hcal_inner, ecal_inner, ecal_outer, hcal_length, ecal_length, radiationLength, interactionLength, inertiatensor, append

using LCIO
using DataStructures

#TODO: procure the correct constants for the ILC
const hcal_outer = 2495
const hcal_inner = 1419
const hcal_length = 6036
const numfaces = 8

const ecal_inner = 1265
const ecal_outer = 1403.5
const ecal_length = 3530
const radiationLength = 25.695807
const interactionLength = 1.041129

"Append contents to the a file with the given filename."
function append(filename, contents)
    file = open(filename, "a")
    write(file, contents)
    close(file)
end

function inertiatensor(positions::Matrix)
    m = zeros(Float64, 3, 3)
    for j=1:size(positions)[2]
        pos = positions[:,j]
        m[1, 1] += (pos[2]^2 + pos[3]^2)*pos[4]
        m[2, 2] += (pos[1]^2 + pos[3]^2)*pos[4]
        m[3, 3] += (pos[1]^2 + pos[2]^2)*pos[4]
        m[1, 2] -= pos[1]*pos[2]*pos[4]
        m[2, 3] -= pos[2]*pos[3]*pos[4]
        m[3, 1] -= pos[1]*pos[3]*pos[4]
    end
    m[2, 1] = m[1, 2]
    m[3, 2] = m[2, 3]
    m[1, 3] = m[3, 1]

    inertia = Symmetric(m)
    return inertia
end

"""
  getPositions(fileName, collectionName, eventNum=1)

Get the position/energy matrix of calorimeter hits in the slcio file, where each output column is [x, y, z, E].
The collection to read from is `collectionName`, and the event index is `eventNum`.
Currently only works for PandoraPFOCollection.
"""
function getpositions(filename, collectionName, eventNum=1)
    for (idx, event) in enumerate(LCIO.open(filename))
        #println("idx: $(idx), eventNum: $(eventNum), equal: $(idx == eventNum)")
        if (idx == eventNum)
            for name in LCIO.getCollectionNameArray(event)
                if (name == collectionName)
                    # collection = LCIO.getCollection(event, name)
                    collection = LCIO.getCollection(event, "PandoraPFOCollection")
                    positions = Array(Float64, 4, 0)

                    for (hitnum, c) in enumerate(collection)
                        for hit in LCIO.getParticleHits(c)
                            positions = hcat(positions, [hit.x, hit.y, hit.z, hit.E])
                        end
                    end
                    return positions
                end
            end
            println("No collection named $(collectionName)")
            throw(KeyError(collectionName))
        end
    end
end

"""
  mapcollections(f, filename, startevent, endevent, collections...)

Map a function f(positionvec4, i) over the specified collections from the specified
start to end event.
"""
function mapcollections(f, filename, startevent, endevent, collections...)
    for (i, event) in zip(1:endevent, LCIO.open(filename))
        positions = nil(Vector{Float64})
        if i < startevent continue end
        for name in LCIO.getCollectionNameArray(event)
            for tname in collections
                if name == tname
                    collection = LCIO.getCollection(event, tname)
                    for hit in collection
                        positions = cons([hit.x, hit.y, hit.z, hit.E], positions)
                    end
                end
            end
        end
        f(hcat(positions...), i)
    end
end

"""
  mappositions(f, filename, startevent, endevent, collectionname="PandoraPFOCollection")

Applies the function `f(positionvec4, i, maxenergy, momentum, particletype)` to the matrix of position-energy 4-vectors,
event number `i`, reconstructed energy `maxenergy`, 3-momentum, and particle type for events `startevent:endevent`.
`positionvec4` is a 4xn matrix where n is the number of hits, and the last component is the energy.
Currently only works for PandoraPFOCollection. Only considers one particle per event (with maximum energy)
"""
function mappositions(f, filename, startevent, endevent, collectionname="PandoraPFOCollection")
    for (i, event) in zip(1:endevent, LCIO.open(filename))
        if i < startevent continue end
        for name in LCIO.getCollectionNameArray(event)
            if name == collectionname
                collection = LCIO.getCollection(event, collectionname)
                positions = nil(Vector{Float64})
                maxEnergy = 0.0
                momentum = Float64[0, 0, 0]
                particletype = 0
                for c in collection
                    p4vec = LCIO.getP4(c)
                    if p4vec.t > maxEnergy
                        particletype = getParticleID(c)
                        maxEnergy = p4vec.t
                        momentum = [p4vec.x, p4vec.y, p4vec.z]
                    end
                    for hit in LCIO.getParticleHits(c)
                        positions = cons([hit.x, hit.y, hit.z, hit.E], positions)
                    end
                end
                len = length(positions)
                positionsvec = Array(Float64, 4, len)
                for (j, pos) in enumerate(positions)
                    positionsvec[:,j] = pos
                end
                f(positionsvec, i, maxEnergy, momentum, particletype)
            end
        end
    end
end

"""
  mapparticles(f, filename, startevent, endevent, collectionname="PandoraPFOCollection")

Applies the function `f(positionvec4, i, maxenergy, momentum, particletype)` to the matrix of position-energy 4-vectors,
event number `i`, reconstructed energy `maxenergy`, momentum, and particle type for all particles in events `startevent:endevent`.
`positionvec4` is a 4xn matrix where n is the number of hits, and the last component is the energy.
Currently only works for PandoraPFOCollection.
"""
function mapparticles(f, filename, startevent, endevent, collectionname="PandoraPFOCollection")
    for (i, event) in zip(1:endevent, LCIO.open(filename))
        if i < startevent continue end
        for name in LCIO.getCollectionNameArray(event)
            if name == collectionname
                collection = LCIO.getCollection(event, collectionname)
                for c in collection
                    positions = nil(Vector{Float64})
                    p4vec = LCIO.getP4(c)
                    maxEnergy = p4vec.t
                    momentum = [p4vec.x, p4vec.y, p4vec.z]
                    particletype = getParticleID(c)
                    for hit in LCIO.getParticleHits(c)
                        positions = cons([hit.x, hit.y, hit.z, hit.E], positions)
                    end
                    len = length(positions)
                    positionsvec = Array(Float64, 4, len)
                    for (j, pos) in enumerate(positions)
                        positionsvec[:,j] = pos
                    end
                    f(positionsvec, i, maxEnergy, momentum, particletype)
                end
            end
        end
    end
end

"""
  maptruthparticles(f, filename, startevent, endevent)

Applies the function `f(positionvec4, i, maxenergy, momentum, particletype)` to the matrix of position-energy 4-vectors,
event number `i`, reconstructed energy `maxenergy`, momentum, and particle type for all particles in events `startevent:endevent`.
Hits are assigned to their true particles, with energies scaled to account for shared hits.
`positionvec4` is a 4xn matrix where n is the number of hits, and the last component is the energy.
"""
function maptruthparticles(f, filename, startevent, endevent)
    for (idx, i) in zip(1:endevent, LCIO.open(filename))
        if idx < startevent continue end
        particleHits = Dict{Ptr{Void}, Dict}()
        ee = getCollection(i, "EcalEndcapHits")
        eb = getCollection(i, "EcalBarrelHits")
        he = getCollection(i, "HcalEndcapHits")
        hb = getCollection(i, "HcalBarrelHits")
        for relation in getCollection(i, "CalorimeterHitRelations")
        	calibratedHit = getRelationFrom(relation)
        	simCaloHit = getRelationTo(relation)
        	pList = getHitMCParticles(simCaloHit)
        	p4 = getCaloHit(calibratedHit)
        	energyList = getHitEnergyList(simCaloHit)
        	simCaloHitEnergy = sum(energyList)
        	# TODO map the particles to the Hits, but correct the hit energy by the sampling fraction
        	# maybe look at energy spread within the cluster
            if (simCaloHit in eb) || (simCaloHit in ee)
                collectionName = "ECAL"
            elseif (simCaloHit in hb) || (simCaloHit in he)
                collectionName = "HCAL"
            else
                # skip muonhits
                continue
            end
            for (kdx, particle) in enumerate(pList)
                h = LCIO.CalHit(p4, energyList[kdx]*p4.E/simCaloHitEnergy)
                if haskey(particleHits, particle)
                    push!(particleHits[particle][collectionName], h)
                else
                    particleHits[particle] = Dict("ECAL" => LCIO.CalHit[], "HCAL" => LCIO.CalHit[])
                end
            end
        end

        for (p, hitDictionary) in particleHits
            p4 = getMCP4(p)
            pdg = getMCPDGid(p)
            len1 = length(hitDictionary["ECAL"])
            len2 = length(hitDictionary["HCAL"])
            positions = Array(Float64, 4, len1+len2)
            ind = 1
            for hit in hitDictionary["ECAL"]
                positions[:, ind] = [hit.x, hit.y, hit.z, hit.E]
                ind += 1
            end
            for hit in hitDictionary["HCAL"]
                positions[:, ind] = [hit.x, hit.y, hit.z, hit.E]
                ind += 1
            end
            f(positions, idx, p4.t, [p4.x, p4.y, p4.z], pdg)
        end
    end

end
end
