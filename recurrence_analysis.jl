using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Random")
Pkg.add("LaTeXStrings")
Pkg.add("Sundials")
Pkg.add("DelimitedFiles")
Pkg.add("Statistics")
Pkg.add("DynamicalSystems")

using DifferentialEquations
using Random
using LaTeXStrings
using Sundials
using DelimitedFiles
using Statistics
using DynamicalSystems

diretory=readdir()

for p in 1:length(diretory)

    #####################Variables S_Max########################

    StatsBlock         = 2     # Microstates size
    rng                = MersenneTwister()
    samples            = 10000 # Samples of Recurrence Microstates
    Parts              = 100   # Threshold
    Low_Eps            = 0.0001
    Var_Eps            = (1.000-Low_Eps)/Parts
    Max_Micro          = Int64(2^(StatsBlock*StatsBlock))
    S_Max              = zeros(Float64,Parts)
    Threshold          = zeros(Float64,Parts)
    Add                = zeros(Int64,Parts)
    a_binary           = zeros(Int64,Parts)
    Stats              = zeros(Float64,Max_Micro*Parts)

    #####################Time and Window########################

    Serie              = readdlm(diretory[p]) #Serie
    t_end,N_of_Comp    = size(Serie)
    Component          = 2
    #t_end             = 170001

    W_Frac             = 0.20 # Fraction of window serie size
    Size               = t_end
    Window_Size        = floor(Int64,W_Frac*Size)
    Jump               = Window_Size # Spaces between windows

    #####################Arguments S_Max########################

    pow_vec           = zeros(Int64,(StatsBlock*StatsBlock))
    for i=1:(StatsBlock*StatsBlock)
        pow_vec[i]=Int64(2^(i-1))
    end

    x_rand            = zeros(Int64,(samples))
    y_rand            = zeros(Int64,(samples))

    for count=1:samples
        x_rand[count]=round(Int64,(rand(rng)*(Window_Size-StatsBlock)))
        y_rand[count]=round(Int64,(rand(rng)*x_rand[count]))
    end

    count_n3            = 0
    count_n2            = 0
    while count_n2 <= (t_end-Window_Size)
        global count_n2=count_n2
        global count_n3=count_n3
        count_n2=count_n2+Jump
        count_n3=count_n3+1
    end

    ############################################################

    function Max_Entropy(Serie)

        Stats[:].=0
        S_Max[:].=0.0

        for i=1:Parts

            Threshold[i]=Low_Eps+(i-1)*Var_Eps

            for count=1:samples
                Add[i]=0
                for count_y=1:StatsBlock
                    for count_x=1:StatsBlock
                        if (abs(Serie[x_rand[count]+count_x]-Serie[y_rand[count]+count_y]) <= Threshold[i])
                            a_binary[i]=1
                        else
                            a_binary[i]=0
                        end

                        Add[i]=Add[i]+a_binary[i]*pow_vec[count_x+((count_y-1)*StatsBlock)]

                    end
                end

                Stats[Add[i]+1+((i-1)*Max_Micro)]+=1

            end

            for j=1:Max_Micro
                if (Stats[j+((i-1)*Max_Micro)] > 0)
                    S_Max[i]+=(-(Stats[j+((i-1)*Max_Micro)]/(1.0*samples))*(log((Stats[j+((i-1)*Max_Micro)]/(1.0*samples)))))
                end
            end

        end

        S_M,Threshold_M=findmax(S_Max)
        Threshold_M=Low_Eps+(Threshold_M-1)*Var_Eps

        return S_M,Threshold_M

    end

    #####################Main Function########################

    max_loops            = count_n3

    V_Mean_Window        = zeros(Float64,Window_Size)

    A1                   = zeros(Float64,max_loops)
    A2                   = zeros(Float64,max_loops)
    A3                   = zeros(Float64,max_loops)
    A4                   = zeros(Float64,max_loops)
    A5                   = zeros(Float64,max_loops)
    A6                   = zeros(Float64,max_loops)
    A7                   = zeros(Float64,max_loops)
    A8                   = zeros(Float64,max_loops)
    A9                   = zeros(Float64,max_loops)
    A10                  = zeros(Float64,max_loops)
    A11                  = zeros(Float64,max_loops)
    A12                  = zeros(Float64,max_loops)
    A13                  = zeros(Float64,max_loops)
    A14                  = zeros(Float64,max_loops)
    S_Vector             = zeros(Float64,max_loops)
    Eps_Vector           = zeros(Float64,max_loops)

    All                  = zeros(Float64,16)
    lmin_v               = 2

    for count2=1:max_loops

        println(count2,' ',max_loops,' ',t_end)

        V_Mean_Window[:]        = Serie[(1+((count2-1)*Jump)):(((count2-1)*Jump)+Window_Size),Component]

        Maximum_Value,Local_Number=findmax(V_Mean_Window)
        Minimum_Value,Local_Number=findmin(V_Mean_Window)
        if ((Maximum_Value-Minimum_Value) != 0.0)
            V_Mean_Window.=((V_Mean_Window.-Minimum_Value)./(Maximum_Value-Minimum_Value))
        end

        S_Vector[count2],Eps_Vector[count2] = Max_Entropy(V_Mean_Window) # For threshold linked to max microstates entropy
        #Eps_Vector[count2]                 = 0.1 #Fixed threshold

        R                                   = RecurrenceMatrix(V_Mean_Window, Eps_Vector[count2]; metric = "euclidean")
        #R                                  = RecurrenceMatrix(V_Mean_Window, 0.1; fixedrate=true, metric = "euclidean") # For threshold linked to recurrence rate

        A1[count2]                          = recurrencerate(R) #recurrence rate

        A2[count2]                          = determinism(R; lmin=lmin_v) #determinism
        A3[count2]                          = dl_average(R; lmin=lmin_v) #average diagonal
        A4[count2]                          = dl_max(R; lmin=lmin_v) # max diagonal
        A5[count2]                          = dl_entropy(R; lmin=lmin_v) # diagonal entropy
        A6[count2]                          = divergence(R) #divergence
        A7[count2]                          = trend(R) #trend

        A8[count2]                          = laminarity(R; lmin=lmin_v) #laminarity
        A9[count2]                          = trappingtime(R; lmin=lmin_v) #trapping time
        A10[count2]                         = vl_average(R; lmin=lmin_v) # average vertical
        A11[count2]                         = vl_max(R; lmin=lmin_v) # max vertical

        A12[count2]                         = meanrecurrencetime(R; lmin=lmin_v) # mean recurrence time
        A13[count2]                         = rt_entropy(R; lmin=lmin_v) #entropy recurrence time
        A14[count2]                         = rt_average(R; lmin=lmin_v) #average mean recurrence time

    end

    if (max_loops < 2)
        All[1]=A1[1];   All[2]=A2[1];   All[3]=A3[1];   All[4]=A4[1];   All[5]=A5[1];
        All[6]=A6[1];   All[7]=A7[1];   All[8]=A8[1];   All[9]=A9[1];   All[10]=A10[1];
        All[11]=A11[1]; All[12]=A12[1]; All[13]=A13[1]; All[14]=A14[1]; All[15]=S_Vector[1];
        All[16]=Eps_Vector[1]
    else
        All[1]=mean(A1);   All[2]=mean(A2);   All[3]=mean(A3);   All[4]=mean(A4);   All[5]=mean(A5);
        All[6]=mean(A6);   All[7]=mean(A7);   All[8]=mean(A8);   All[9]=mean(A9);   All[10]=mean(A10);
        All[11]=mean(A11); All[12]=mean(A12); All[13]=mean(A13); All[14]=mean(A14); All[15]=mean(S_Vector);
        All[16]=mean(Eps_Vector)

        All_Std  = zeros(Float64,16)
        All_Std[1]=std(A1);   All_Std[2]=std(A2);   All_Std[3]=std(A3);   All_Std[4]=std(A4);   All_Std[5]=std(A5);
        All_Std[6]=std(A6);   All_Std[7]=std(A7);   All_Std[8]=std(A8);   All_Std[9]=std(A9);   All_Std[10]=std(A10);
        All_Std[11]=std(A11); All_Std[12]=std(A12); All_Std[13]=std(A13); All_Std[14]=std(A14); All_Std[15]=std(S_Vector);
        All_Std[16]=std(Eps_Vector)
        b=string("z_Out__Std_",a[p])
        writedlm(b, All_Std)
    end

    c=string("z_Out_", a[p])

    writedlm(c, All)
end