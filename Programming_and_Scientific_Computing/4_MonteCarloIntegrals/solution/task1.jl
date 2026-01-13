using Random
using Statistics
using Makie, GLMakie

# Compute the monte carlo weights based on the function f(x) = x * sin(x)
# Note that, π is the normalization factor for uniform sampling in [0, π]
function f(x)
    return π * x * sin(x)
end

function sample_x(n_samples::Int)
    samples = rand(n_samples) * π  # Sample uniformly in [0, π]
    return samples
end

function monte_carlo_integration(f, n_samples::Int)
    samples = sample_x(n_samples)

    # Monte carlo weights
    ωi = f.(samples)
    avgω = mean(ωi)

    # Compute standard error
    s2 = sum((ωi .- avgω) .^2) / (n_samples - 1)

    return avgω, sqrt(s2/n_samples)
end

function plot_convergence(f, samples, exact_value)
    estimates = Float64[]
    errors = Float64[]

    for n in samples
        estimate, error = monte_carlo_integration(f, n)
        push!(estimates, estimate)
        push!(errors, 1.96*error)
    end

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Number of Samples", ylabel = "Integral Estimate", title = "Monte Carlo Integration Convergence", xscale = log10)
    xlims!(ax, samples[1], samples[end])
    scatter!(ax, samples, estimates, label = "Estimate", color = :blue)
    band!(ax, samples, estimates .- errors, estimates .+ errors, color = (:blue, 0.3), label = "Error Band")
    hlines!(ax, [exact_value], color = :red, linestyle = :dash, label = "Exact Value")
    axislegend(ax)
    fig
end

function plot_hist(f, N, samples)
    estimates = Float64[]
    errors = Float64[]

    for _ in 1:samples
        estimate, error = monte_carlo_integration(f, N)
        push!(estimates, estimate)
        push!(errors, error)
    end

    fig = Figure(resolution = (800, 600))

    # Plot results as a histogram
    ax = Axis(fig[1, 1], xlabel = "Integral Estimate", ylabel = "Frequency", title = "Histogram of Monte Carlo Estimates")
    hist!(ax, estimates, bins = 30, color = :green, label = "Estimates")

    std_from_errors = mean(errors)

    # Print average and standard deviation
    avg_estimate = mean(estimates) 
    std_estimate = std(estimates)
    vlines!(ax, [avg_estimate], color = :red, linestyle = :dash, label = "Average Estimate")
    # Lines for std deviation
    vlines!(ax, [avg_estimate - std_estimate, avg_estimate + std_estimate], color = :orange, linestyle = :dot, label = "Std Deviation")
    vlines!(ax, [avg_estimate - std_from_errors, avg_estimate + std_from_errors], color = :purple, linestyle = :dashdot, label = "Avg Error Std Dev")
    axislegend(ax)
    fig 
end

function error_vs_N(f, N_values, n_trials, exact_value)
    avg_errors = Float64[]

    for N in N_values
        errors = Float64[]
        for _ in 1:n_trials
            estimate, error = monte_carlo_integration(f, N)
            push!(errors, abs(estimate - exact_value))
        end
        push!(avg_errors, mean(errors))
    end

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Number of Samples", ylabel = "Average Absolute Error", title = "Error vs Number of Samples", xscale = log10, yscale = log10)
    scatter!(ax, N_values, avg_errors, color = :magenta, label = "Avg Absolute Error")

    # Plot linear fit
    log_N = log10.(N_values)
    log_errors = log10.(avg_errors) 

    X = hcat(ones(length(log_N)), log_N)
    coeffs = X \ log_errors
    fit_errors = 10 .^ (X * coeffs)

    lines!(ax, N_values, fit_errors, color = :black, label = "Linear Fit")

    text!(ax, "Fit: Error ≈ N^$(round(coeffs[2], digits=2))", position = (0.2, 0.8), space=:relative, color = :black)
    axislegend(ax)
    fig
end