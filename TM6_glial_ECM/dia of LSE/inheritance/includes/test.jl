function check_inheritance_vertical(u0_start, u0_end, len)

    for index in range(1, len - 1, step = 1)

        if u0_end[1, index, :] != u0_start[1, index+1, :]
            println("Test inheritance vertical FAIELD");
            println("Index: $(index)");
            return -1;
        end
    end

    println("Test inheritance vertical SUCCEEDED");
    return 0;
end
