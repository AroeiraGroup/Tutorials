using Random
using Statistics    
using Makie, GLMakie

function uniform_monte_carlo_weights(N)
    samples = rand(N) * π  # Sample uniformly in [0, π]
    ωi = [π * x * sin(x) for x in samples]
    return ωi
end

function importance_sampling_weights(N)
    U = rand(N)  # Sample uniformly in [0, 1]
    X = acos.(1 .- 2 .* U)  # Inverse transform sampling for p(x) = (2/π) * sin^2(x)
    return 2*X
end

function plot_comparison(Nvalues; exact_value=π)

    uniform_values = Float64[]
    uniform_errors = Float64[]
    importance_values = Float64[]
    importance_errors = Float64[]

    for n in Nvalues
        # Uniform sampling
        ωi_uniform = uniform_monte_carlo_weights(n)
        avg_uniform = mean(ωi_uniform)
        s2_uniform = sum((ωi_uniform .- avg_uniform) .^2) / (n - 1)
        push!(uniform_values, avg_uniform)
        push!(uniform_errors, sqrt(s2_uniform/n))

        # Importance sampling
        ωi_importance = importance_sampling_weights(n)
        avg_importance = mean(ωi_importance)
        s2_importance = sum((ωi_importance .- avg_importance) .^2) / (n - 1)
        push!(importance_values, avg_importance)
        push!(importance_errors, sqrt(s2_importance/n))
    end

    # Plot results
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Number of Samples", ylabel = "Integral Estimate", title = "Uniform vs Importance Sampling", xscale = log10)
    xlims!(ax, Nvalues[1], Nvalues[end])   
    scatter!(ax, Nvalues, uniform_values, label = "Uniform Sampling", color = :blue)
    scatter!(ax, Nvalues, importance_values, label = "Importance Sampling", color = :green)
    hlines!(ax, [exact_value], color = :red, linestyle = :dash, label = "Exact Value")

    # Plot errors as bands
    band!(ax, Nvalues, uniform_values .- 1.96 .* uniform_errors, uniform_values .+ 1.96 .* uniform_errors, color = (:blue, 0.3), label = "Uniform 95% CI")
    band!(ax, Nvalues, importance_values .- 1.96 .* importance_errors, importance_values .+ 1.96 .* importance_errors, color = (:green, 0.3), label = "Importance 95% CI")
    axislegend(ax)
    fig
end

function plot_variance_ratios(Nvalues)
    uniform_variances = Float64[]
    importance_variances = Float64[]

    for n in Nvalues
        # Uniform sampling
        ωi_uniform = uniform_monte_carlo_weights(n)
        s2_uniform = sum((ωi_uniform .- mean(ωi_uniform)) .^2) / (n - 1)
        push!(uniform_variances, s2_uniform)

        # Importance sampling
        ωi_importance = importance_sampling_weights(n)
        s2_importance = sum((ωi_importance .- mean(ωi_importance)) .^2) / (n - 1)
        push!(importance_variances, s2_importance)
    end

    variance_ratios = importance_variances ./ uniform_variances

    # Plot variance ratios
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Number of Samples", ylabel = "Variance Ratio (Importance / Uniform)", title = "Variance Ratio vs Number of Samples", xscale = log10)
    xlims!(ax, Nvalues[1], Nvalues[end])
    scatter!(ax, Nvalues, variance_ratios, label = "Variance Ratio", color = :purple)
    axislegend(ax)
    fig
end