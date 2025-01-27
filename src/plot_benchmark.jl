@recipe function f(benchmark::BenchmarkResults; title = "", color = :auto, label = "", x = :callcount, showribbon = true)

    markersize := 3
    titlefontsize := 10
  
    @series begin
        seriestype := :line
        color := color
        label := label
        subplot := 1
        xlabel := x == :callcount ? "Function calls" : "Iterations"
        ylabel := "Success rate"
        ylim := (0,1)
        
        title := title
    
        if x == :callcount
            x = benchmark.callcount
        elseif x == :run_length 
            x = benchmark.run_length
        else
            error("Valid options for x are : :run_length or :callcount")
        end
        y = benchmark.success_rate

        yerror = (abs.(benchmark.success_rate_qlow .- y), benchmark.success_rate_qhigh .- y)
        if showribbon
            ribbon := yerror
        end

        x, y
    end

end