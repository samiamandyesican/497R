using Plots

function generate_plots(labels, data, foldername)
    j = 1
    mkpath("data/$foldername/V_th")
    mkpath("data/$foldername/V_tv")

    for i in labels
        
        if i == "V_th"
            j += 1
            println("skip success")

        elseif i == "V_tv"
            j += 1
            println("skip success")

        else

            rm("data/$foldername/V_th/Vth_v_$(labels[j]).png", force=true)
            rm("data/$foldername/V_tv/Vtv_v_$(labels[j]).png", force=true)

            plt1 = plot(data[:, 1], data[:, j], xlabel=labels[1], ylabel=labels[j])
            savefig(plt1, "data/$foldername/V_th/Vth_v_$(labels[j]).png")
            println("plot $(labels[1]) vs $(labels[j]) success")

            plt2 = plot(data[:, 2], data[:, j], xlabel=labels[2], ylabel=labels[j])
            savefig(plt2, "data/$foldername/V_tv/Vtv_v_$(labels[j]).png")
            println("plot $(labels[2]) vs $(labels[j]) success")

            j +=1

        end

    end

end