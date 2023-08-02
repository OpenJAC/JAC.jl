
    struct CascadeData{T}
        lines::Array{T,1}
    end
    
    wa = CascadeData{Int64}([1,2,3])
    wb = CascadeData{Float64}([1.,2.,3.])
    
    wc = CascadeData[]
    push!(wc, wa);    push!(wc, wb);  
    
    for data in wc
        @show data.lines[1]
    end
    
    @show wc[1].lines[1]
    
