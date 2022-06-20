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
        # println(lineSplit)
        chromRegionPosition = string(lineSplit[1]) * "_" * string(lineSplit[2]) * "_" * string(lineSplit[3])
        # println(bed2Data[chromRegionPosition])
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

#@code_native f(upperLim)
#@code_lowered f(upperLim)  # I love this, very similar to python bytecode
#@time     # also wonderful for timing pieces of code


# def parse_beds():
#     bedlist1 = []
#     bedlist2 = []
#     bed1count = 0
#     bed2count = 0
#     matches_conditions = 0
#     with open(bed1,'r') as inF1:
#         lines = inF1.readlines()
#         for line in lines:
#             line = line.split()
#             chrom1, start1, end1, TE_name = str(line[0]),int(line[1]),int(line[2]),str(line[3])
#             bed1count += 1
#             with open(bed2,'r') as inF2:
#                 linez = inF2.readlines()
#                 for lin in linez:
#                     lin = lin.split()
#                     chrom2, start2, end2, score2 = str(lin[0]),int(lin[1]),int(lin[2]),float(lin[3])
#                     bed2count += 1
#                     #assumes a perfect match between file 1 & 2 for chromosome,position and region
#                     if chrom1 == chrom2 and start1 and start2 >= start1 and end1 and end2 <= end1:
#                         #the append is an ordered append, adding the scores inline with one another for downstream correlations
#                        # bedlist1.append(score1)
#                         #bedlist2.append(score2)
#                         matches_conditions += 1
#                         print(chrom1, start1, end1, TE_name, score2, file=open("resultsfile.txt", "a"))

#     combined_list = [bedlist1,bedlist2]
#     print("[STDOUT]: Example Bedlists:")
#     print("[STDOUT]: BedList1:",bedlist1[0:5],"\n[STDOUT]: Bedlist2:",bedlist2[0:5])
#     print("[STDOUT]: Out of:",bed1count,"regions in --bed1",float(100*(matches_conditions/bed1count)),"% or",matches_conditions,"exact matching regions were found in --bed2")
#     print("[STDOUT]: The correlation of the measured values in these regions are as follows:")
#     return combined_list