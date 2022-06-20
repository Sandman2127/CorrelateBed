using ArgParse
using Statistics

function parse_commandline()
    s = ArgParseSettings(description = "Correlate the region values in column 5 of a bed file against another bedfile. i.e. H3K36me3 vs DNA methylation")
    @add_arg_table! s begin
        "--bed1"
        help = "A normal bedfile you want to correlate regions against bed2, requires columns: chromosome\tstart\tend\tanyValue\tscore"
        "--bed2"
        help = "A normal bedfile you want to correlate regions against bed1, requires columns: chromosome\tstart\tend\tanyValue\tscore"
    end
    return parse_args(s)
end

function parseBed2(bedPath)
    lines = readlines(bedPath)
    Bed2Data = Dict{String}{Float64}()
    for line in lines
        lineSplit = split(line)
        key = string(lineSplit[1]) * "_" * string(lineSplit[2]) * "_" * string(lineSplit[3])
        #println(lineSplit[5])
        Bed2Data[key] = parse(Float64,lineSplit[5])
    end
    return Bed2Data
end

function findMatches(bed1Path,bed2Data,foutName)
    bed1CorrVals = Float64[]
    bed2CorrVals = Float64[]
    count = 0
    foundMatch = 0
    for line in eachline(bed1Path)
        count += 1 
        lineSplit = split(line)
        chromRegionPosition = string(lineSplit[1]) * "_" * string(lineSplit[2]) * "_" * string(lineSplit[3])
        if haskey(bed2Data,chromRegionPosition)
            foundMatch += 1
            push!(bed1CorrVals,parse(Float64,lineSplit[5]))
            push!(bed2CorrVals,bed2Data[chromRegionPosition])
        end
    end
    # provide clues as to which data was correlated and how much of it matched
    percentFound = 100*Float64(foundMatch/count)
    println("[STDOUT]: out of the total of ",count," lines in bed1, we saw perfect matches at:",foundMatch," positions or a total of: ",percentFound,"%")
    println("[STDOUT]: printing first 5 values of bed1 correlations")
    println(bed1CorrVals[1:5])
    println("[STDOUT]: printing first 5 values of bed2 correlations")
    println(bed2CorrVals[1:5])
    output = Any[bed1CorrVals,bed2CorrVals]
    return output
end

function main()
    parsed_args = parse_commandline()
    bed2Data = parseBed2(parsed_args["bed2"])
    foutName = parsed_args["bed1"][1:end-4] * ".intersect.bed"
    corVals = findMatches(parsed_args["bed1"],bed2Data,foutName)
    bed1CorVals,bed2CorVals = corVals[1],corVals[2]
    println("[STDOUT]: the correlation between bed1:",parsed_args["bed1"]," and bed2: ",parsed_args["bed2"]," is: ",Statistics.cor(bed1CorVals,bed2CorVals))
end

@time main()