#' Plotting for REMIND-DIETER coupling
#'
#' Read in REMIND and DIETER results directly from gdx files and create REMIND-DIETER_validation.pdf
#'
#' @param outputdir path to the output directory of the REMIND-DIETER coupled run
#'
#' @author Adrian Odenweller
#' @import dplyr
#' @import ggplot2
#' @importFrom stringr str_sort
#' @importFrom lusweave swopen swlatex swclose
#' @importFrom quitte read.gdx revalue.levels
#' @importFrom gridExtra arrangeGrob
#' @importFrom tidyr gather spread
#' @export
plotDIETERvalidation <- function(outputdir) {
    # Configurations ----------------------------------------------------------

    report.periods <- c(seq(2015, 2060, 5), seq(2070, 2100, 10))

    remind.nonvre.mapping <- c(
        coalchp = "Coal (Lig + HC)",
        igcc = "Coal (Lig + HC)",
        igccc = "Coal (Lig + HC)",
        pcc = "Coal (Lig + HC)",
        pco = "Coal (Lig + HC)",
        pc = "Coal (Lig + HC)",
        tnrs = "Nuclear",
        ngt = "OCGT",
        ngcc = "CCGT",
        ngccc = "CCGT",
        gaschp = "CCGT",
        biochp = "Biomass",
        bioigcc = "Biomass",
        bioigccc = "Biomass",
        NULL
    )

    remind.vre.mapping <- c(
        hydro = "Hydro",
        wind = "Wind",
        spv = "Solar",
        NULL
    )

    remind.tech.mapping <- c(remind.nonvre.mapping, remind.vre.mapping)

    dieter.tech.exclude <- c("OCGT_ineff", "Wind_off")

    dieter.tech.mapping <- c(
        hc = "Hard coal",
        lig = "Lignite",
        coal = "Coal (Lig + HC)",
        nuc = "Nuclear",
        OCGT_eff = "OCGT",
        CCGT = "CCGT",
        bio = "Biomass",
        ror = "Hydro",
        Wind_on = "Wind",
        Solar = "Solar",
        NULL
    )

    color.mapping <- c(
        "CCGT" = "#999959", "Lignite" = "#0c0c0c", "Coal (Lig + HC)" = "#0c0c0c",
        "Solar" = "#ffcc00", "Wind" = "#337fff", "Biomass" = "#005900",
        "OCGT" = "#e51900", "Hydro" = "#191999", "Nuclear" = "#ff33ff",
        "Hard coal" = "#808080"
    )

    sm_TWa_2_MWh <- 8.76E9

    # Directories -------------------------------------------------------------

    report.output.file <- file.path(outputdir, "REMIND-DIETER_validation.pdf")

    remind.files <- list.files(outputdir, pattern = "fulldata_[0-9]+\\.gdx") %>%
        str_sort(numeric = TRUE)
    cat(paste0("REMIND files: ", length(remind.files), "\n"))

    dieter.files <- list.files(outputdir, pattern = "results_DIETER_i[0-9]+\\.gdx") %>%
        str_sort(numeric = TRUE)
    cat(paste0("DIETER files: ", length(dieter.files), "\n"))

    dieter.files.report <- list.files(outputdir, pattern = "report_DIETER_i[0-9]+\\.gdx") %>%
        str_sort(numeric = TRUE)
    cat(paste0("DIETER report files: ", length(dieter.files.report), "\n"))

    # Determine iteration step of DIETER
    dieter.iter.step <- floor(length(remind.files) / length(dieter.files))
    cat(paste0("DIETER iter step: ", dieter.iter.step, "\n"))

    # Define functions ----------------------------------------------------

    DIETERplotCapacityFactors <- function() {
        all_regi <- all_te <- cap <- capfac <- char <- NULL
        generation <- iteration <- model <- rlf <- NULL
        ttot <- variable <- value <- tall <- vm_cap <- vm_capFac <- technology <- var <- NULL

        # Data preparation (REMIND) -----------------------------------------------

        cat("Plot capacity factors \n")

        out.remind.capfac <- NULL
        for (i in 1:length(remind.files)) {

            # Capacity factor for non-VRE ---------------------------------------------

            # Read in vm_capFac
            remind.vm_capFac <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("vm_capFac", fields = "l", squeeze = F) %>%
                rename(tall = ttot) %>%
                mutate(variable = "vm_capFac")

            # Read in vm_cap
            remind.vm_cap <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("vm_cap", fields = "l", squeeze = F) %>%
                filter(rlf == 1) %>%
                select(-rlf) %>%
                mutate(variable = "vm_cap")

            # Join both vm_capFac and vm_cap
            remind.data.nonVRE <- rbind(remind.vm_capFac, remind.vm_cap) %>%
                spread(variable, value) %>%
                filter(all_regi == "DEU") %>%
                filter(tall %in% report.periods) %>%
                filter(all_te %in% names(remind.nonvre.mapping)) %>%
                mutate(technology = all_te) %>%
                revalue.levels(technology = remind.nonvre.mapping) %>%
                mutate(generation = vm_cap * vm_capFac) %>%
                group_by(tall, technology) %>%
                summarise(capfac = sum(generation) / sum(vm_cap)) %>%
                mutate(iteration = i)

            # Capacity factor for VRE -------------------------------------------------

            # Read in pm_dataren with VRE capacity factors over grades
            remind.pm_dataren <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("pm_dataren", squeeze = F) %>%
                filter(all_te %in% names(remind.vre.mapping)) %>%
                filter(char == "nur") %>%
                select(-char) %>%
                rename(capfac = value)

            # Read in vm_capDistr with VRE capacity distribution over grades
            remind.vm_capDistr <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("vm_capDistr", fields = "l", squeeze = F) %>%
                rename(cap = value)

            # Join pm_dataren with vm_capDistr and calculate VRE CFs
            remind.data.VRE <- left_join(remind.pm_dataren, remind.vm_capDistr) %>%
                filter(all_regi == "DEU") %>%
                filter(tall %in% report.periods) %>%
                rename(technology = all_te) %>%
                revalue.levels(technology = remind.vre.mapping) %>%
                mutate(generation = cap * capfac) %>%
                group_by(tall, technology) %>%
                summarise(capfac = sum(generation) / sum(cap)) %>%
                mutate(iteration = i)

            out.remind.capfac <- rbind(out.remind.capfac, remind.data.nonVRE, remind.data.VRE) %>%
                mutate(model = "REMIND")
        }

        # Data preparation (DIETER) -----------------------------------------------

        out.dieter.capfac <- NULL
        for (i in 1:length(dieter.files)) {
            dieter.data <- file.path(outputdir, dieter.files[i]) %>%
                read.gdx("report4RM", squeeze = F, colNames = c("file", "tall", "all_regi", "technology", "var", "value")) %>%
                select(!c(file, all_regi)) %>%
                filter(tall %in% report.periods) %>%
                mutate(tall = as.numeric(as.character(tall))) %>%
                filter(var == "capfac") %>%
                filter(!technology %in% dieter.tech.exclude) %>%
                revalue.levels(technology = dieter.tech.mapping) %>%
                mutate(iteration = dieter.iter.step * i) %>%
                mutate(model = "DIETER")

            out.dieter.capfac <- rbind(out.dieter.capfac, dieter.data)
        }

        # Make capacity factors available in parent environment
        out.dieter.capfac <<- out.dieter.capfac
        out.remind.capfac <<- out.remind.capfac

        # Plotting ----------------------------------------------------------------

        swlatex(sw, paste0("\\section{Capacity factors}"))

        for (t.rep in report.periods) {
            plot.remind <- out.remind.capfac %>%
                filter(tall == t.rep)

            plot.dieter <- out.dieter.capfac %>%
                filter(tall == t.rep) %>%
                filter(!technology %in% c("Lignite", "Hard coal"))

            swlatex(sw, paste0("\\subsection{Capacity factors in ", t.rep, "}"))

            p <- ggplot() +
                geom_line(data = plot.remind, aes(x = iteration, y = capfac, color = model)) +
                geom_point(data = plot.dieter, aes(x = iteration, y = value, color = model)) +
                xlab("Iteration") +
                ylab("Capacity factor") +
                facet_wrap(~technology, nrow = 3)

            swfigure(sw, print, p)
        }

        swlatex(sw, "\\subsection{Capacity factors over time (last iteration)}")

        plot.remind <- out.remind.capfac %>%
            filter(iteration == max(iteration))

        plot.dieter <- out.dieter.capfac %>%
            filter(iteration == max(iteration)) %>%
            filter(!technology %in% c("Lignite", "Hard coal"))

        p <- ggplot() +
            geom_line(data = plot.remind, aes(x = tall, y = capfac, color = model)) +
            geom_line(data = plot.dieter, aes(x = tall, y = value, color = model)) +
            facet_wrap(~technology, nrow = 3) +
            xlab("Time") +
            ylab("Capacity factor")

        swfigure(sw, print, p)
    }

    DIETERplotCapacities <- function() {
        all_regi <- all_te <- capacity <- demand <- iteration <- NULL
        rlf <- tall <- technology <- value <- var <- NULL

        # Data preparation (REMIND) -----------------------------------------------

        cat("Plot capacities \n")

        out.remind.cap <- NULL
        out.remind.dem <- NULL
        for (i in 1:length(remind.files)) {

            # Read in vm_cap (capacity)
            remind.vm_cap <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("vm_cap", field = "l", squeeze = F) %>%
                filter(rlf == 1) %>%
                select(-rlf) %>%
                filter(all_regi == "DEU") %>%
                filter(tall %in% report.periods) %>%
                filter(all_te %in% names(remind.tech.mapping)) %>%
                revalue.levels(all_te = remind.tech.mapping) %>%
                mutate(all_te = factor(all_te, levels = rev(unique(remind.tech.mapping)))) %>%
                group_by(tall, all_te) %>%
                summarise(capacity = 1e3 * sum(value)) %>% # REMIND capacity is in TW
                ungroup() %>%
                mutate(iteration = i)

            # Read in v32_seelDem (total secondary electricity demand)
            remind.v32_seelDem <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("v32_seelDem", field = "l", squeeze = F) %>%
                filter(all_regi == "DEU") %>%
                filter(tall %in% report.periods) %>%
                mutate(demand = 0.000155891 * 8760 * 1e3 * value) %>%
                mutate(iteration = i)

            out.remind.cap <- rbind(out.remind.cap, remind.vm_cap)
            out.remind.dem <- rbind(out.remind.dem, remind.v32_seelDem)
        }


        # Data preparation (DIETER) -----------------------------------------------

        out.dieter <- NULL
        for (i in 1:length(dieter.files)) {
            dieter.data <- file.path(outputdir, dieter.files[i]) %>%
                read.gdx("report4RM", squeeze = F, colNames = c("file", "tall", "all_regi", "technology", "var", "value")) %>%
                select(tall, technology, var, value) %>%
                filter(tall %in% report.periods) %>%
                mutate(tall = as.numeric(as.character(tall))) %>%
                filter(var == "capacity") %>%
                mutate(capacity = value / 1e3) %>% # DIETER capacity is in MW
                select(-var) %>%
                filter(!technology %in% dieter.tech.exclude) %>%
                revalue.levels(technology = dieter.tech.mapping) %>%
                mutate(technology = factor(technology, levels = rev(unique(dieter.tech.mapping)))) %>%
                mutate(iteration = dieter.iter.step * i)

            out.dieter <- rbind(out.dieter, dieter.data)
        }

        # Plotting ----------------------------------------------------------------

        swlatex(sw, paste0("\\section{Capacities}"))

        for (t.rep in report.periods) {
            plot.remind.cap <- out.remind.cap %>%
                filter(tall == t.rep)

            plot.remind.dem <- out.remind.dem %>%
                filter(tall == t.rep)

            plot.dieter <- out.dieter %>%
                filter(tall == t.rep)

            swlatex(sw, paste0("\\subsection{Capacities in ", t.rep, "}"))

            p1 <- ggplot() +
                geom_area(data = plot.remind.cap, aes(x = iteration, y = capacity, fill = all_te), alpha = 0.5) +
                scale_fill_manual(name = "Technology", values = color.mapping) +
                geom_line(data = plot.remind.dem, aes(x = iteration, y = demand, colour = "Peak demand"), linetype = "dotted") +
                scale_color_manual(name = "Demand", values = c("Peak demand" = "black")) +
                coord_cartesian(xlim = c(0, max(plot.remind.cap$iteration))) +
                xlab("Iteration") +
                ylab("Capacity [GW]") +
                ggtitle("REMIND")

            p2 <- ggplot() +
                geom_area(data = plot.dieter, aes(x = iteration, y = capacity, fill = technology), alpha = 0.5) +
                scale_fill_manual(name = "Technology", values = color.mapping) +
                geom_line(data = plot.remind.dem, aes(x = iteration, y = demand, colour = "Peak demand"), linetype = "dotted") +
                scale_color_manual(name = "Demand", values = c("Peak demand" = "black")) +
                coord_cartesian(xlim = c(0, max(plot.remind.cap$iteration))) +
                xlab("Iteration") +
                ylab("Capacity [GW]") +
                ggtitle("DIETER")

            grid.newpage()
            p <- arrangeGrob(rbind(ggplot2::ggplotGrob(p1), ggplot2::ggplotGrob(p2)))

            swfigure(sw, grid.draw, p)
        }


        swlatex(sw, "\\subsection{Capacities over time (last iteration)}")

        plot.remind.cap <- out.remind.cap %>%
            filter(iteration == max(out.remind.cap$iteration))

        plot.remind.dem <- out.remind.dem %>%
            filter(iteration == max(out.remind.dem$iteration))

        plot.dieter <- out.dieter %>%
            filter(iteration == max(out.dieter$iteration))

        p1 <- ggplot() +
            geom_area(data = plot.remind.cap, aes(x = tall, y = capacity, fill = all_te), alpha = 0.5) +
            scale_fill_manual(name = "Technology", values = color.mapping) +
            geom_line(data = plot.remind.dem, aes(x = tall, y = demand, colour = "Peak demand"), linetype = "dotted") +
            scale_color_manual(name = "Demand", values = c("Peak demand" = "black")) +
            xlab("Time") +
            ylab("Capacity [GW]") +
            ggtitle("REMIND")

        p2 <- ggplot() +
            geom_area(data = plot.dieter, aes(x = tall, y = capacity, fill = technology), alpha = 0.5) +
            scale_fill_manual(name = "Technology", values = color.mapping) +
            geom_line(data = plot.remind.dem, aes(x = tall, y = demand, colour = "Peak demand"), linetype = "dotted") +
            scale_color_manual(name = "Demand", values = c("Peak demand" = "black")) +
            xlab("Time") +
            ylab("Capacity [GW]") +
            ggtitle("DIETER")

        grid.newpage()
        p <- arrangeGrob(rbind(ggplot2::ggplotGrob(p1), ggplot2::ggplotGrob(p2)))

        swfigure(sw, grid.draw, p)
    }

    DIETERplotGeneration <- function() {
        all_enty.1 <- all_regi <- tall <- all_te <- value <- NULL
        var <- technology <- iteration <- generation <- NULL

        # Data preparation (REMIND) -----------------------------------------------

        cat("Plot generation \n")

        out.remind <- NULL
        for (i in 1:length(remind.files)) {

            # Read in vm_prodSe (generation)
            remind.vm_prodSe <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("vm_prodSe", field = "l", squeeze = F) %>%
                filter(all_enty.1 == "seel") %>%
                filter(all_regi == "DEU") %>%
                filter(tall %in% report.periods) %>%
                filter(all_te %in% names(remind.tech.mapping)) %>%
                revalue.levels(all_te = remind.tech.mapping) %>%
                mutate(all_te = factor(all_te, levels = rev(unique(remind.tech.mapping)))) %>%
                group_by(tall, all_te) %>%
                summarise(generation = 8760 * sum(value)) %>% # TWa ->TWh
                ungroup() %>%
                mutate(iteration = i)

            out.remind <- rbind(out.remind, remind.vm_prodSe)
        }

        # Data preparation (DIETER) -----------------------------------------------

        out.dieter <- NULL
        for (i in 1:length(dieter.files)) {
            dieter.data <- file.path(outputdir, dieter.files[i]) %>%
                read.gdx("report4RM", squeeze = F, colNames = c("file", "tall", "all_regi", "technology", "var", "value")) %>%
                select(tall, technology, var, value) %>%
                mutate(tall = as.numeric(as.character(tall))) %>%
                filter(tall %in% report.periods) %>%
                filter(var == "generation") %>%
                mutate(generation = value / 1e6) %>% # MWh -> TWh
                select(-var) %>%
                filter(!technology %in% dieter.tech.exclude) %>%
                revalue.levels(technology = dieter.tech.mapping) %>%
                mutate(technology = factor(technology, levels = rev(unique(dieter.tech.mapping)))) %>%
                mutate(iteration = dieter.iter.step * i)

            out.dieter <- rbind(out.dieter, dieter.data)
        }

        # Plotting ----------------------------------------------------------------

        swlatex(sw, paste0("\\section{Generation}"))

        for (t.rep in report.periods) {
            plot.remind <- out.remind %>%
                filter(tall == t.rep)

            plot.dieter <- out.dieter %>%
                filter(tall == t.rep)

            swlatex(sw, paste0("\\subsection{Generation in ", t.rep, "}"))

            p1 <- ggplot() +
                geom_area(data = plot.remind, aes(x = iteration, y = generation, fill = all_te), alpha = 0.5) +
                scale_fill_manual(name = "Technology", values = color.mapping) +
                coord_cartesian(xlim = c(0, max(plot.remind$iteration))) +
                xlab("Iteration") +
                ylab("Generation [TWh]") +
                ggtitle("REMIND")

            p2 <- ggplot() +
                geom_area(data = plot.dieter, aes(x = iteration, y = generation, fill = technology), alpha = 0.5) +
                scale_fill_manual(name = "Technology", values = color.mapping) +
                coord_cartesian(xlim = c(0, max(plot.remind$iteration))) +
                xlab("Iteration") +
                ylab("Generation [TWh]") +
                ggtitle("DIETER")

            grid.newpage()
            p <- arrangeGrob(rbind(ggplot2::ggplotGrob(p1), ggplot2::ggplotGrob(p2)))

            swfigure(sw, grid.draw, p)
        }


        swlatex(sw, "\\subsection{Generation over time (last iteration)}")

        plot.remind <- out.remind %>%
            filter(iteration == max(out.remind$iteration))

        plot.dieter <- out.dieter %>%
            filter(iteration == max(out.dieter$iteration))

        p1 <- ggplot() +
            geom_area(data = plot.remind, aes(x = tall, y = generation, fill = all_te), alpha = 0.5) +
            scale_fill_manual(name = "Technology", values = color.mapping) +
            xlab("Time") +
            ylab("Generation [TWh]") +
            ggtitle("REMIND")

        p2 <- ggplot() +
            geom_area(data = plot.dieter, aes(x = tall, y = generation, fill = technology), alpha = 0.5) +
            scale_fill_manual(name = "Technology", values = color.mapping) +
            xlab("Time") +
            ylab("Generation [TWh]") +
            ggtitle("DIETER")

        grid.newpage()
        p <- arrangeGrob(rbind(ggplot2::ggplotGrob(p1), ggplot2::ggplotGrob(p2)))

        swfigure(sw, grid.draw, p)
    }

    DIETERplotAddedCapacities <- function() {
        var <- technology <- variable <- model <- iteration <- NULL
        tall <- value <- NULL
        # Data preparation --------------------------------------------------------

        cat("Plot added capacities \n")

        dieter.report.cap <- c("DIETER pre-investment capacities", "REMIND pre-investment capacities")
        dieter.report.addcap <- c("DIETER added capacities (GW)", "REMIND added capacities (GW)")
        dieter.report.divest <- c("REMIND divestment (GW)")

        dieter.report.vars <- c(dieter.report.cap, dieter.report.addcap, dieter.report.divest)

        out.dieter.report <- NULL
        for (i in 1:length(dieter.files.report)) {
            dieter.data <- file.path(outputdir, dieter.files.report[i]) %>%
                read.gdx("report_tech", squeeze = F, colNames = c("file", "model", "tall", "all_regi", "var", "technology", "value")) %>%
                filter(var %in% dieter.report.vars) %>%
                filter(!technology == "coal") %>%
                revalue.levels(technology = dieter.tech.mapping) %>%
                mutate(variable = case_when(
                    var %in% dieter.report.cap ~ "Pre-inv. cap.",
                    var %in% dieter.report.addcap ~ "Added cap.",
                    var %in% dieter.report.divest ~ "Divestment"
                )) %>%
                mutate(variable = factor(variable, levels = rev(c("Pre-inv. cap.", "Added cap.", "Divestment")))) %>%
                mutate(iteration = dieter.iter.step * i)

            out.dieter.report <- rbind(out.dieter.report, dieter.data)
        }


        # Plotting ----------------------------------------------------------------

        swlatex(sw, "\\onecolumn")
        swlatex(sw, paste0("\\section{Added capacities}"))

        for (iter.rep in 1:length(dieter.files.report)) {
            swlatex(sw, paste0("\\subsection{Added capacities in iteration ", iter.rep * dieter.iter.step, "}"))

            plot.dieter <- out.dieter.report %>%
                filter(model == "DIETER") %>%
                filter(iteration == iter.rep * dieter.iter.step) %>%
                mutate(tall = as.numeric(as.character(tall)) - 1) # Shift for dodged plot

            plot.remind <- out.dieter.report %>%
                filter(model == "REMIND") %>%
                filter(iteration == iter.rep * dieter.iter.step) %>%
                mutate(tall = as.numeric(as.character(tall)) + 1) %>% # Shift for dodged plot
                mutate(value = ifelse(variable == "Divestment", -value, value)) # Divestment has negative value

            p <- ggplot() +
                geom_bar(data = plot.dieter, aes(x = tall, y = value, fill = model, alpha = variable), colour = "black", stat = "identity", position = "stack", width = 2) +
                geom_bar(data = plot.remind, aes(x = tall, y = value, fill = model, alpha = variable), colour = "black", stat = "identity", position = "stack", width = 2) +
                scale_alpha_manual(values = c("Pre-inv. cap." = 1, "Added cap." = 0.5, "Divestment" = 0.2), limits = c("Pre-inv. cap.", "Added cap.", "Divestment")) +
                facet_wrap(~technology, scales = "free") +
                coord_cartesian(xlim = c(2010, 2100)) +
                theme(legend.position = "bottom") +
                xlab("Time") +
                ylab("Capacity [GW]")

            swfigure(sw, print, p, sw_option = "width=20, height=10")
        }
        swlatex(sw, "\\twocolumn")
    }

    DIETERplotLCOEs <- function(out.dieter.capfac) {
        var <- tall <- technology <- iteration <- var <- cost <- model <- value <- NULL
        IC <- OM <- FC <- CO2 <- NULL
        `CO2 cost` <- `O&M cost` <- `annualized investment cost` <- NULL
        capfac <- `fuel cost (divided by eta)` <- files <- all_regi <- NULL

        # Data preparation --------------------------------------------------------

        cat("Plot LCOEs \n")

        dieter.report.lcoe.kW <- c("annualized investment cost", "O&M cost")
        dieter.report.lcoe.MWh <- c("fuel cost (divided by eta)", "CO2 cost")
        dieter.report.mv <- "DIETER Market value ($/MWh)"

        dieter.report.vars <- c(dieter.report.lcoe.kW, dieter.report.lcoe.MWh, dieter.report.mv)

        out.dieter.report <- NULL
        for (i in 1:length(dieter.files.report)) {
            dieter.data <- file.path(outputdir, dieter.files.report[i]) %>%
                read.gdx("report_tech", squeeze = F, colNames = c("file", "model", "tall", "all_regi", "var", "technology", "value")) %>%
                filter(var %in% dieter.report.vars) %>%
                filter(tall %in% report.periods) %>%
                mutate(var = factor(var, levels = rev(dieter.report.vars))) %>%
                filter(!technology %in% c(dieter.tech.exclude, "coal")) %>%
                revalue.levels(technology = dieter.tech.mapping) %>%
                mutate(iteration = dieter.iter.step * i)

            out.dieter.report <- rbind(out.dieter.report, dieter.data)
        }

        out.dieter.lcoe <- out.dieter.report %>%
            filter(var %in% c(dieter.report.lcoe.kW, dieter.report.lcoe.MWh)) %>%
            select(!c(file, all_regi)) %>%
            rbind(out.dieter.capfac) %>%
            spread(var, value) %>%
            mutate(IC = 1e3 * `annualized investment cost` / (capfac * 8760)) %>%
            mutate(OM = 1e3 * `O&M cost` / (capfac * 8760)) %>%
            mutate(FC = `fuel cost (divided by eta)`) %>%
            mutate(CO2 = `CO2 cost`) %>%
            gather(c("IC", "OM", "FC", "CO2"), key = "var", value = "cost") %>%
            select(c(model, tall, technology, iteration, var, cost)) %>%
            filter(!technology == "Coal (Lig + HC)") %>%
            mutate(var = factor(var, levels = rev(c("IC", "OM", "FC", "CO2"))))

        out.dieter.mv <- out.dieter.report %>%
            filter(var %in% dieter.report.mv)


        # Plotting ----------------------------------------------------------------

        swlatex(sw, "\\onecolumn")
        swlatex(sw, paste0("\\section{Levelised cost of electricity (LCOE)}"))

        for (t.rep in report.periods) {
            swlatex(sw, paste0("\\subsection{LCOEs in ", t.rep, "}"))

            plot.dieter.lcoe <- out.dieter.lcoe %>%
                filter(tall == t.rep)

            plot.dieter.mv <- out.dieter.mv %>%
                filter(tall == t.rep)

            plot.dieter.capfac <- out.dieter.capfac %>%
                filter(tall == t.rep) %>%
                filter(!technology == "Coal (Lig + HC)")

            p <- ggplot() +
                geom_bar(data = plot.dieter.lcoe, aes(x = iteration, y = cost, fill = var), stat = "identity", position = "stack") +
                scale_fill_discrete(name = "LCOE components", labels = c("Annualised investment", "O&M", "Fuel", "CO2"), limits = c("IC", "OM", "FC", "CO2")) +
                geom_point(data = plot.dieter.mv, aes(x = iteration, y = value, colour = "Market value")) +
                scale_color_manual(values = "black", name = NULL) +
                geom_text(data = plot.dieter.capfac, aes(x = iteration, y = 350, label = paste0(100 * round(value, 2), "%"))) +
                theme(legend.position = "bottom") +
                xlab("Iteration") +
                ylab("LCOE [$/MWh]") +
                coord_cartesian(ylim = c(0, 400)) +
                facet_wrap(~technology)

            swfigure(sw, print, p, sw_option = "width=20, height=10")
        }

        swlatex(sw, paste0("\\subsection{LCOEs over time (last iteration)}"))

        plot.dieter.lcoe <- out.dieter.lcoe %>%
            filter(iteration == max(iteration)) %>%
            mutate(tall = as.numeric(as.character(tall)))

        plot.dieter.mv <- out.dieter.mv %>%
            filter(iteration == max(iteration)) %>%
            mutate(tall = as.numeric(as.character(tall)))

        plot.dieter.capfac <- out.dieter.capfac %>%
            filter(iteration == max(iteration)) %>%
            filter(tall %in% report.periods) %>%
            filter(!technology == "Coal (Lig + HC)")

        p <- ggplot() +
            geom_bar(data = plot.dieter.lcoe, aes(x = tall, y = cost, fill = var), stat = "identity", position = "stack") +
            scale_fill_discrete(name = "LCOE components", labels = c("Annualised investment", "O&M", "Fuel", "CO2"), limits = c("IC", "OM", "FC", "CO2")) +
            geom_point(data = plot.dieter.mv, aes(x = tall, y = value, colour = "Market value")) +
            scale_color_manual(values = "black", name = NULL) +
            geom_text(data = plot.dieter.capfac, aes(x = tall, y = 350, label = paste0(100 * round(value, 2), "%"))) +
            theme(legend.position = "bottom") +
            xlab("Time") +
            ylab("LCOE [$/MWh]") +
            coord_cartesian(ylim = c(0, 400)) +
            facet_wrap(~technology, scale = "free")

        swfigure(sw, print, p, sw_option = "width=20, height=10")

        swlatex(sw, "\\twocolumn")
    }

    DIETERplotSeelPrice <- function(out.remind.capfac) {
        all_regi <- ttot <- all_enty <- tall <- q32_balSe.m <- m <- NULL
        qm_budget.m <- seel.price <- iteration <- capfac <- technology <- NULL


        # Data preparation --------------------------------------------------------

        cat("Plot Seel price \n")

        out.remind.seel <- NULL
        for (i in 1:length(remind.files)) {
            remind.q32_balSe <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("q32_balSe", field = "m", squeeze = F) %>%
                filter(all_regi == "DEU") %>%
                select(!all_regi) %>%
                filter(ttot %in% report.periods) %>%
                select(!all_enty) %>%
                rename(tall = ttot) %>%
                rename(q32_balSe.m = m)

            remind.qm_budget <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("qm_budget", field = "m", squeeze = F) %>%
                filter(all_regi == "DEU") %>%
                select(!all_regi) %>%
                filter(ttot %in% report.periods) %>%
                rename(tall = ttot) %>%
                rename(qm_budget.m = m)

            remind.seel <- left_join(remind.q32_balSe, remind.qm_budget) %>%
                mutate(seel.price = 1.2 * 1e12 / sm_TWa_2_MWh * q32_balSe.m / qm_budget.m) %>% # (10^12 2005$)/TWa -> 2015$/MWh
                mutate(iteration = i)

            out.remind.seel <- rbind(out.remind.seel, remind.seel)
        }


        # Plotting ----------------------------------------------------------------

        swlatex(sw, "\\onecolumn")
        swlatex(sw, paste0("\\section{Secondary electricity (seel) price + capacity factors}"))

        swlatex(sw, paste0("\\subsection{Seel price and capacity factor over iterations}"))

        p <- ggplot() +
            geom_line(data = out.remind.seel, aes(x = iteration, y = seel.price, size = "Seel")) +
            scale_size_manual(name = "Shadow price", values = 1) +
            geom_line(data = out.remind.capfac, aes(x = iteration, y = 2.5 * 100 * capfac, colour = technology)) +
            scale_colour_manual(name = "Capacity factor", values = color.mapping) +
            scale_y_continuous(name = "Seel price [$/MWh]", limits = c(0, 250), sec.axis = sec_axis(~ . / 2.5, name = paste0("CF", "(%)"))) +
            theme(legend.position = "bottom") +
            xlab("Iteration") +
            facet_wrap(~tall, ncol = 4)

        swfigure(sw, print, p, sw_option = "width=20, height=12")

        swlatex(sw, "\\twocolumn")
        swlatex(sw, paste0("\\subsection{Seel price and capacity factor over time (last iteration)}"))

        plot.remind.seel <- out.remind.seel %>%
            filter(iteration == max(iteration))

        plot.remind.capfac <- out.remind.capfac %>%
            filter(iteration == max(iteration))

        p <- ggplot() +
            geom_line(data = plot.remind.seel, aes(x = tall, y = seel.price, size = "Seel")) +
            scale_size_manual(name = "Shadow price", values = 1) +
            geom_line(data = plot.remind.capfac, aes(x = tall, y = 2.5 * 100 * capfac, colour = technology)) +
            scale_colour_manual(name = "Capacity factor", values = color.mapping) +
            scale_y_continuous(name = "Seel price [$/MWh]", limits = c(0, 250), sec.axis = sec_axis(~ . / 2.5, name = paste0("CF", "(%)"))) +
            theme(legend.position = "bottom") +
            xlab("Time")

        swfigure(sw, print, p)
    }

    DIETERplotPeakDemandPrice <- function() {

        all_enty <- tall <- q32_peakDemand_DT.m <- m <- NULL
        all_regi <- ttot <- tall <- qm_budget.m <- NULL
        peakdem.price <- iteration <- NULL

        # Data preparation --------------------------------------------------------

        cat("Plot peak demand price \n")

        out.remind.peakdem <- NULL
        for (i in 1:length(remind.files)) {
            remind.q32_peakDemand_DT <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("q32_peakDemand_DT", field = "m", squeeze = F) %>%
                select(!all_enty) %>%
                filter(tall %in% report.periods) %>%
                rename(q32_peakDemand_DT.m = m)

            remind.qm_budget <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("qm_budget", field = "m", squeeze = F) %>%
                filter(all_regi == "DEU") %>%
                select(!all_regi) %>%
                filter(ttot %in% report.periods) %>%
                rename(tall = ttot) %>%
                rename(qm_budget.m = m)

            remind.peakdem <- left_join(remind.q32_peakDemand_DT, remind.qm_budget) %>%
                mutate(peakdem.price = 1.2 * 1e3 * q32_peakDemand_DT.m / qm_budget.m) %>% # (10^12 2005$)/TW-> 2015$/kW
                mutate(iteration = i)

            out.remind.peakdem <- rbind(out.remind.peakdem, remind.peakdem)
        }


        # Plotting ----------------------------------------------------------------

        swlatex(sw, "\\onecolumn")
        swlatex(sw, paste0("\\section{Peak demand (capital constraint) shadow price}"))

        swlatex(sw, paste0("\\subsection{Peak demand shadow price over iterations}"))

        p <- ggplot() +
            geom_line(data = out.remind.peakdem, aes(x = iteration, y = peakdem.price, size = "Peak demand (capital constraint)")) +
            scale_size_manual(name = "Shadow price", values = 1) +
            theme(legend.position = "bottom") +
            xlab("Iteration") +
            ylab("Shadow price [$/kW]") +
            facet_wrap(~tall, ncol = 4)

        swfigure(sw, print, p, sw_option = "width=20, height=12")

        swlatex(sw, "\\twocolumn")
        swlatex(sw, paste0("\\subsection{Peak demand shadow price over time (last iteration)}"))

        plot.remind.peakdem <- out.remind.peakdem %>%
            filter(iteration == max(iteration))

        p <- ggplot() +
            geom_line(data = plot.remind.peakdem, aes(x = tall, y = peakdem.price, size = "Peak demand (capital constraint)")) +
            scale_size_manual(name = "Shadow price", values = 1) +
            theme(legend.position = "bottom") +
            xlab("Time") +
            ylab("Shadow price [$/kW]")

        p

        swfigure(sw, print, p)
    }

    # LaTeX configurations ----------------------------------------------------

    template <- c(
        "\\documentclass[a4paper,landscape,twocolumn]{article}",
        "\\setlength{\\oddsidemargin}{-0.8in}",
        "\\setlength{\\evensidemargin}{-0.5in}",
        "\\setlength{\\topmargin}{-0.8in}",
        "\\setlength{\\parindent}{0in}",
        "\\setlength{\\headheight}{0in}",
        "\\setlength{\\topskip}{0in}",
        "\\setlength{\\headsep}{0in}",
        "\\setlength{\\footskip}{0.2in}",
        "\\setlength\\textheight{0.95\\paperheight}",
        "\\setlength\\textwidth{0.95\\paperwidth}",
        "\\setlength{\\parindent}{0in}",
        "\\usepackage{float}",
        "\\usepackage[bookmarksopenlevel=section,colorlinks=true,linkbordercolor={0.9882353 0.8352941 0.7098039}]{hyperref}",
        "\\hypersetup{bookmarks=true,pdfauthor={PIK}}",
        "\\usepackage{graphicx}",
        "\\usepackage[strings]{underscore}",
        "\\usepackage{Sweave}",
        "\\begin{document}",
        "<<echo=false>>=",
        "options(width=110)",
        "@"
    )

    # Open LaTeX PDF ----------------------------------------------------------

    sw <- swopen(report.output.file, template = template)

    swlatex(sw, "\\tableofcontents\\newpage")

    # Capacity factors --------------------------------------------------------

    DIETERplotCapacityFactors()

    # Capacities --------------------------------------------------------------

    DIETERplotCapacities()

    # Generation --------------------------------------------------------------

    DIETERplotGeneration()

    # Added capacities --------------------------------------------------------

    DIETERplotAddedCapacities()

    # LCOEs -------------------------------------------------------------------

    DIETERplotLCOEs(out.dieter.capfac)

    # Price: Secondary electricity price --------------------------------------

    DIETERplotSeelPrice(out.remind.capfac)

    # Price: Peak demand ------------------------------------------------------

    DIETERplotPeakDemandPrice()

    # (Residual) load duration curves -----------------------------------------

    # DIETERplotRLDCs()

    # Price duration curves ---------------------------------------------------

    # DIETERplotPriceDurationCurves()

    # (Inverse) screening curves ----------------------------------------------

    # DIETERplotInverseScreeningCurve()

    # Markups -----------------------------------------------------------------


    # Close LaTeX PDF ---------------------------------------------------------

    swclose(sw)
}