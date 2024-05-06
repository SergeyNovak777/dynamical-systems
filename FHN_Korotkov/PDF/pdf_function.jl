function get_peaks(data; level_zero = "above")
    xrange, trange = data
    if level_zero == "above"
        xrange, trange = xrange[findall(xrange.>=0)], trange[findall(xrange.>=0)]
    elseif level_zero == "below"
        xrange, trange = xrange[findall(xrange.<=0)], trange[findall(xrange.<=0)]
    elseif  level_zero =="nothing"
        nothing
    end
    array_peaks = Float64[]
    array_t_peaks = Float64[]
    for index in range(2, length(xrange)-1 )
        if xrange[index-1] < xrange[index] > xrange[index+1]
            push!(array_peaks, xrange[index])
            push!(array_t_peaks, trange[index])
        end
    end
    return array_peaks, array_t_peaks
end

function get_peaks_neg(data)
    xrange, trange = data
    xrange, trange = xrange[findall(xrange.<=0)], trange[findall(xrange.<=0)]
    array_peaks = Float64[]
    array_t_peaks = Float64[]
    for index in range(2, length(xrange)-1 )
        if xrange[index-1] > xrange[index] < xrange[index+1]
            push!(array_peaks, xrange[index])
            push!(array_t_peaks, trange[index])
        end
    end
    return array_peaks, array_t_peaks
end

Hs_above(x, k) = Statistics.mean(x) + k * Statistics.std(x)
Hs_below(x, k) = Statistics.mean(x) - k * Statistics.std(x)

function pdf(array_peaks, count_thesholds, maxvalue, ϵ)
    threshold_range = range(0.0, maxvalue, length = count_thesholds)
    total_number_peaks = length(array_peaks)
    array_PDF = Float64[]
    for index in range(1, count_thesholds, step = 1)
        number_peaks_crossing_theshold = length(findall( (threshold_range[index] + ϵ) .>= array_peaks .>= threshold_range[index] ) )
        PDF_for_theshold = number_peaks_crossing_theshold / total_number_peaks
        push!(array_PDF, PDF_for_theshold)
    end
    return threshold_range, array_PDF
end

function pdf_v1(array_peaks, count_thesholds, maxvalue, ϵ)
    threshold_range = range(0.0, maxvalue, length = count_thesholds)
    total_number_peaks = length(array_peaks)
    array_PDF = similar(array_peaks, Float64, count_thesholds)  # Preallocate array_PDF

    threshold_plus_epsilon = threshold_range .+ ϵ  # Calculate threshold + ϵ outside the loop

    for index in 1:count_thesholds
        number_peaks_crossing_theshold = count((threshold_plus_epsilon[index] .>= array_peaks) .& (array_peaks .>= threshold_range[index]))
        PDF_for_theshold = number_peaks_crossing_theshold / total_number_peaks
        array_PDF[index] = PDF_for_theshold
    end

    return threshold_range, array_PDF
end


function pdf_v2(array_peaks, count_thesholds, minvalue, maxvalue, ϵ)
    threshold_range = range(minvalue, maxvalue, length = count_thesholds)
    total_number_peaks = length(array_peaks)
    array_PDF = similar(array_peaks, Float64, count_thesholds)  # Preallocate array_PDF

    threshold_plus_epsilon = threshold_range .+ ϵ  # Calculate threshold + ϵ outside the loop

    for index in 1:count_thesholds
        number_peaks_crossing_theshold = count((threshold_plus_epsilon[index] .>= array_peaks) .& (array_peaks .>= threshold_range[index]))
        PDF_for_theshold = number_peaks_crossing_theshold / total_number_peaks
        array_PDF[index] = PDF_for_theshold
    end

    return threshold_range, array_PDF
end

# старый вариант
function CALCPDF(spikes, count_thesholds)
    thesholds = range(0,maximum(spikes),count_thesholds)
    ee_counter = [sum(i-> s<=i<s+0.9, spikes) for s in thesholds]
    pdf = ee_counter ./ length(spikes)
    return thesholds, pdf
end