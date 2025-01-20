@recipe function f(benchmark::BenchmarkResults; title = nothing, color = :auto, label = "")

    markersize := 3
    titlefontsize := 10
  
    @series begin
        seriestype := :line
        color := color
        label := label
        subplot := 1
        xlabel := "Run length"
        ylabel := "Success rate"
        ylim := (0,1)
        if title != nothing
            title := title
        end

        x = benchmark.run_length
        y = benchmark.success_rate

        yerror = (abs.(benchmark.success_rate_qlow .- y), benchmark.success_rate_qhigh .- y)
        ribbon := yerror

        x, y
    end

end