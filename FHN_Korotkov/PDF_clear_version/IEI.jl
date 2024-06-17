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

function get_PDF(events, shift)
    count_events = length(events)
    PDFs = zeros(count_events)

    for index in range(1, count_events, step = 1)
        count_event_i = count(events[index]-shift .<= events .<= events[index]+shift)
        PDF_event_i = count_event_i / count_events
        PDFs[index] = PDF_event_i
    end
    return PDFs
end

function get_PDF_without_less_shift(events, shift)
    count_events = length(events)
    PDFs = zeros(count_events)

    for index in range(1, count_events, step = 1)
        count_event_i = count(events[index].<= events .<= events[index]+shift)
        PDF_event_i = count_event_i / count_events
        PDFs[index] = PDF_event_i
    end
    return PDFs
end

function CALCPDF(spikes, thesholds)
    #thesholds = range(minimum(spikes),maximum(spikes),count_thesholds)
    ee_counter = [sum(i-> s<=i<s+025, spikes) for s in thesholds]
    pdf = ee_counter ./ length(spikes)
    return pdf
end