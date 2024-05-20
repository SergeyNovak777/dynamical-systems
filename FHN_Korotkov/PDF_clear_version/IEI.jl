function get_IEI(t_peaks)
    count_IEI = length(t_peaks)-1
    array_IEI = Float64[]
    for index in 1:count_IEI
        t_peak_i = t_peaks[index]
        t_peak_next = t_peaks[index+1]
        IEI = abs(t_peak_next - t_peak_i)
        push!(array_IEI, IEI)
    end
    return array_IEI
end

function get_PDF_IEI(IEI; shift = 10)
    total_count_IEI = length(IEI)
    array_PDF = Float64[]
    for index in 1:total_count_IEI
        count_IEI_i = count(IEI[index]-shift .<= IEI .<= IEI[index]+shift)
        PDF_IEI_i = count_IEI_i / total_count_IEI
        push!(array_PDF, PDF_IEI_i)
    end
    return array_PDF
end