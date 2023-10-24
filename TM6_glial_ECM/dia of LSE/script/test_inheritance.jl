function get_log_inheritance_vertical(u0s, u0ttrs, len)
    file_logs = "log_inheritance_vertical $(len)x$(len).txt";
        for index_int in range(1, len, step = 1)
            open(file_logs, "a") do file
                print(file, "index_int: $(index_int) \n");
                print(file, "u0 start: $(u0ttrs[1, index_int, :]) \n");
                print(file, "u0 end: $(u0s[1, index_int, :]) \n");
                print(file, "\n");
            end
        end
end

function get_log_inheritance_horizontal(u0s, u0ttrs, len)
    file_logs = "log_inheritance_horizontal $(len)x$(len).txt"
    for index_ext in range(1, len, step = 1)
        for index_int in range(1, len, step = 1)
            open(file_logs, "a") do file
                print(file, "index_int: $(index_int) \n");
                print(file, "index_ext: $(index_ext) \n");
                print(file, "u0 start: $(u0ttrs[index_int, index_ext, :]) \n");
                print(file, "u0 end: $(u0s[index_int, index_ext, :]) \n");
                print(file, "\n");
            end
        end
    end 
end

function check_inheritance_vertical(u0s, u0ttrs, len)
    retcode = 0;

    for index in range(1, len - 1, step = 1)
        if u0s[1,index,:]!=u0ttrs[1,index+1,:]
            println("Test inheritance vertical FAIELD");
            println("Index: $(index)");
            retcode = -1;
            break;
        end    
    end

    if retcode == 0
        println("Test inheritance vertical SUCCEEDED");
    end
end

function check_inheritance_horizontal(u0s, u0ttrs, len)
    retcode = 0;
    for index_ext in range(1, len, step = 1)
        for index_int in range(2, len, step = 1)
            if u0s[index_int-1, index_ext, :] == u0ttrs[index_int, index_ext, :]
                nothing
            else
                println("int: $(index_int-1); ext: $(index_ext)");
                println("int: $(index_int); ext: $(index_ext)");
                println("");
                retcode = -1
            end
        end
    end 
    if retcode == 0
        println("Test inheritance horizontal SUCCEEDED");
    end
end

function checks(u0s, u0ttrs)
    len, _, _ = size(u0s);

    #get_log_inheritance_vertical(u0s, u0ttrs, len);
    #get_log_inheritance_horizontal(u0s, u0ttrs, len);

    check_inheritance_vertical(u0s, u0ttrs, len);
    check_inheritance_horizontal(u0s, u0ttrs, len);
end
