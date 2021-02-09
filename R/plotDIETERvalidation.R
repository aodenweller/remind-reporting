#' Plotting for REMIND-DIETER coupling
#'
#' Read in REMIND and DIETER results directly from gdx files and create REMIND-DIETER_validation.pdf
#'
#' @param outputdir path to the output directory of the REMIND-DIETER coupled run
#'
#' @author Adrian Odenweller
#' @importFrom dplyr %>% mutate select filter rename summarise group_by ungroup case_when left_join inner_join full_join desc
#' @importFrom ggplot2 ggplot geom_line aes xlab ylab facet_wrap geom_area scale_fill_manual scale_color_manual coord_cartesian ggtitle geom_bar scale_alpha_manual theme scale_x_continuous scale_y_continuous scale_size_manual guides scale_linetype_manual sec_axis guide_legend
#' @importFrom stringr str_sort
#' @importFrom lusweave swopen swlatex swfigure swclose
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

        # Set variables to NULL for code check compliance
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

        return(list("DIETER" = out.dieter.capfac, "REMIND" = out.remind.capfac))
    }

    DIETERplotCapacities <- function() {

        # Set variables to NULL for code check compliance
        all_regi <- all_te <- capacity <- demand <- iteration <- NULL
        rlf <- tall <- technology <- value <- var <- NULL

        # Data preparation (REMIND) -----------------------------------------------

        cat("Plot capacities \n")

        out.remind.cap <- NULL
        out.remind.dem <- NULL
        for (i in 1:length(remind.files)) {

            # Read in vm_cap (capacity)
            remind.vm_cap <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("vm_cap", fields = "l", squeeze = F) %>%
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
                read.gdx("v32_seelDem", fields = "l", squeeze = F) %>%
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
                geom_line(data = plot.remind.dem, aes(x = iteration, y = demand, color = "Peak demand"), linetype = "dotted") +
                scale_color_manual(name = "Demand", values = c("Peak demand" = "black")) +
                coord_cartesian(xlim = c(0, max(plot.remind.cap$iteration))) +
                xlab("Iteration") +
                ylab("Capacity [GW]") +
                ggtitle("REMIND")

            p2 <- ggplot() +
                geom_area(data = plot.dieter, aes(x = iteration, y = capacity, fill = technology), alpha = 0.5) +
                scale_fill_manual(name = "Technology", values = color.mapping) +
                geom_line(data = plot.remind.dem, aes(x = iteration, y = demand, color = "Peak demand"), linetype = "dotted") +
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
            geom_line(data = plot.remind.dem, aes(x = tall, y = demand, color = "Peak demand"), linetype = "dotted") +
            scale_color_manual(name = "Demand", values = c("Peak demand" = "black")) +
            xlab("Time") +
            ylab("Capacity [GW]") +
            ggtitle("REMIND")

        p2 <- ggplot() +
            geom_area(data = plot.dieter, aes(x = tall, y = capacity, fill = technology), alpha = 0.5) +
            scale_fill_manual(name = "Technology", values = color.mapping) +
            geom_line(data = plot.remind.dem, aes(x = tall, y = demand, color = "Peak demand"), linetype = "dotted") +
            scale_color_manual(name = "Demand", values = c("Peak demand" = "black")) +
            xlab("Time") +
            ylab("Capacity [GW]") +
            ggtitle("DIETER")

        grid.newpage()
        p <- arrangeGrob(rbind(ggplot2::ggplotGrob(p1), ggplot2::ggplotGrob(p2)))

        swfigure(sw, grid.draw, p)
    }

    DIETERplotGeneration <- function() {

        # Set variables to NULL for code check compliance
        all_enty.1 <- all_regi <- tall <- all_te <- value <- NULL
        var <- technology <- iteration <- generation <- NULL

        # Data preparation (REMIND) -----------------------------------------------

        cat("Plot generation \n")

        out.remind <- NULL
        for (i in 1:length(remind.files)) {

            # Read in vm_prodSe (generation)
            remind.vm_prodSe <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("vm_prodSe", fields = "l", squeeze = F) %>%
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

        # Set variables to NULL for code check compliance
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
                read.gdx("report_tech",
                    squeeze = F,
                    colNames = c("file", "model", "tall", "all_regi", "var", "technology", "value")
                ) %>%
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
                geom_bar(data = plot.dieter, aes(x = tall, y = value, fill = model, alpha = variable), color = "black", stat = "identity", position = "stack", width = 2) +
                geom_bar(data = plot.remind, aes(x = tall, y = value, fill = model, alpha = variable), color = "black", stat = "identity", position = "stack", width = 2) +
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

        # Set variables to NULL for code check compliance
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
                read.gdx("report_tech",
                    squeeze = F,
                    colNames = c("file", "model", "tall", "all_regi", "var", "technology", "value")
                ) %>%
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
                geom_point(data = plot.dieter.mv, aes(x = iteration, y = value, color = "Market value")) +
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
            geom_point(data = plot.dieter.mv, aes(x = tall, y = value, color = "Market value")) +
            scale_color_manual(values = "black", name = NULL) +
            geom_text(data = plot.dieter.capfac, aes(x = tall, y = 350, label = paste0(100 * round(value, 2), "%"))) +
            theme(legend.position = "bottom") +
            xlab("Time") +
            ylab("LCOE [$/MWh]") +
            coord_cartesian(ylim = c(0, 400)) +
            facet_wrap(~technology, scales = "free")

        swfigure(sw, print, p, sw_option = "width=20, height=10")

        swlatex(sw, "\\twocolumn")
    }

    DIETERplotSeelPrice <- function(out.remind.capfac) {

        # Set variables to NULL for code check compliance
        all_regi <- ttot <- all_enty <- tall <- q32_balSe.m <- m <- NULL
        qm_budget.m <- seel.price <- iteration <- capfac <- technology <- NULL

        # Data preparation --------------------------------------------------------

        cat("Plot Seel price \n")

        out.remind.seel <- NULL
        for (i in 1:length(remind.files)) {
            remind.q32_balSe <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("q32_balSe", fields = "m", squeeze = F) %>%
                filter(all_regi == "DEU") %>%
                select(!all_regi) %>%
                filter(ttot %in% report.periods) %>%
                select(!all_enty) %>%
                rename(tall = ttot) %>%
                rename(q32_balSe.m = m)

            remind.qm_budget <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("qm_budget", fields = "m", squeeze = F) %>%
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
            geom_line(data = out.remind.capfac, aes(x = iteration, y = 2.5 * 100 * capfac, color = technology)) +
            scale_color_manual(name = "Capacity factor", values = color.mapping) +
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
            geom_line(data = plot.remind.capfac, aes(x = tall, y = 2.5 * 100 * capfac, color = technology)) +
            scale_color_manual(name = "Capacity factor", values = color.mapping) +
            scale_y_continuous(name = "Seel price [$/MWh]", limits = c(0, 250), sec.axis = sec_axis(~ . / 2.5, name = paste0("CF", "(%)"))) +
            theme(legend.position = "bottom") +
            xlab("Time")

        swfigure(sw, print, p)
    }

    DIETERplotPeakDemandPrice <- function() {

        # Set variables to NULL for code check compliance
        all_enty <- tall <- q32_peakDemand_DT.m <- m <- NULL
        all_regi <- ttot <- tall <- qm_budget.m <- NULL
        peakdem.price <- iteration <- NULL

        # Data preparation --------------------------------------------------------

        cat("Plot peak demand price \n")

        out.remind.peakdem <- NULL
        for (i in 1:length(remind.files)) {
            remind.q32_peakDemand_DT <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("q32_peakDemand_DT", fields = "m", squeeze = F) %>%
                select(!all_enty) %>%
                filter(tall %in% report.periods) %>%
                rename(q32_peakDemand_DT.m = m)

            remind.qm_budget <- file.path(outputdir, remind.files[i]) %>%
                read.gdx("qm_budget", fields = "m", squeeze = F) %>%
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

    DIETERplotRLDCs <- function() {

        # Set variables to NULL for code check compliance
        model <- all_regi <- var <- value <- tall <- demand.RLDC <- NULL
        technology <- Solar <- Solar_curtailed <- Solar.RLDC <- Wind <- NULL
        Wind_curtailed <- Wind.RLDC <- hour.sorted <- iteration <- demand.hour.sorted <- NULL


        # Data preparation --------------------------------------------------------

        cat("Plot RLDCs \n")

        # Order of technologies in RLDC plot
        rldc.order <-
            c(
                "Nuclear",
                "Hydro",
                "Lignite",
                "Hard coal",
                "Biomass",
                "CCGT",
                "OCGT"
            )
        vre.order <- c("Wind", "Solar")
        order <- c(rldc.order, vre.order)

        # Initialise output files
        out.dieter.demand <- NULL
        out.dieter.rldc <- NULL
        # Loop over DIETER iterations
        for (i in 1:length(dieter.files.report)) {
            # Read in demand(hour)
            dieter.report_hours <-
                file.path(outputdir, dieter.files.report[i]) %>%
                read.gdx("report_hours",
                    squeeze = F,
                    colNames = c("file", "model", "tall", "all_regi", "var", "hour", "value")
                ) %>%
                select(!c(file, model, all_regi)) %>%
                filter(var == "fixed demand") %>%
                select(!var) %>%
                rename(demand.RLDC = value) %>% # Rename demand for RLDC calculation later
                mutate(hour = as.numeric(substring(hour, 2))) %>%
                group_by(tall) %>%
                arrange(desc(demand.RLDC)) %>% # Descending order
                mutate(demand.hour.sorted = seq(1, 8760)) %>% # Descending hour (sorted)
                mutate(iteration = dieter.iter.step * i) %>%
                ungroup()

            out.dieter.demand <- rbind(out.dieter.demand, dieter.report_hours)

            # Read in generation(hour,tech)
            dieter.report_tech_hours <-
                file.path(outputdir, dieter.files.report[i]) %>%
                read.gdx("report_tech_hours",
                    squeeze = F,
                    colNames = c("file", "model", "tall", "all_regi", "var", "technology", "hour", "value")
                ) %>%
                select(!c(file, model, all_regi)) %>%
                revalue.levels(technology = dieter.tech.mapping) %>%
                mutate(hour = as.numeric(substring(hour, 2))) %>%
                mutate(technology = as.character(technology)) %>%
                mutate(technology = case_when( # Make curtailment a separate "technology"
                    var == "curtailment of fluct res" ~ paste0(technology, "_curtailed"),
                    TRUE ~ technology
                )) %>%
                select(!var) %>%
                complete(technology, tall, hour = 1:8760, fill = list(value = 0)) # Fill up missing hours with 0

            # Join both datasets for RLDC calculation
            dieter.data <-
                inner_join(dieter.report_hours, dieter.report_tech_hours) %>%
                spread(technology, value)

            # Calculate RLDCs for solar and wind (with curtailment)
            dieter.data <- dieter.data %>%
                mutate(Solar.RLDC = demand.RLDC - Solar - Solar_curtailed) %>%
                group_by(tall) %>%
                arrange(desc(Solar.RLDC), group_by = T) %>%
                mutate(Solar.hour.sorted = seq(1, 8760)) %>%
                mutate(Wind.RLDC = Solar.RLDC - Wind - Wind_curtailed) %>%
                arrange(desc(Wind.RLDC), group_by = T) %>%
                mutate(Wind.hour.sorted = seq(1, 8760))

            # Calculate RLDCs for dispatchable technologies (without curtailment)
            # This loop calculates the RLDC lines, with different x-axes for each technology
            vars <- c("Wind", rev(rldc.order))
            for (t in 1:(length(vars) - 1)) {
                var1 <- vars[t]
                var2 <- vars[t + 1]
                dieter.data <- dieter.data %>%
                    mutate(!!paste0(var2, ".RLDC") := !!sym(paste0(var1, ".RLDC")) - !!sym(var2)) %>%
                    arrange(desc(!!sym(paste0(var2, ".RLDC")))) %>%
                    mutate(!!paste0(var2, ".hour.sorted") := seq(1, 8760))
            }

            # Calculate differences between RLDC lines for area plot
            # This loop calculates the height of each stacked sorted technology for the same x-axis
            vars <- c("demand", rev(order))
            for (v in 1:(length(vars) - 1)) {
                var1 <- vars[v]
                var2 <- vars[v + 1]
                for (t in sort(unique(dieter.data$tall))) {
                    dieter.temp <- dieter.data %>%
                        select(
                            tall,
                            paste0(var1, ".RLDC"),
                            paste0(var1, ".hour.sorted"),
                            paste0(var2, ".RLDC"),
                            paste0(var2, ".hour.sorted")
                        ) %>%
                        filter(tall == t)

                    # Sort both RLDCs from large to small for subsequent subtraction
                    m1 <- match(1:8760, dieter.temp[paste0(var1, ".hour.sorted")][[1]])
                    m2 <- match(1:8760, dieter.temp[paste0(var2, ".hour.sorted")][[1]])
                    rldc1 <- dieter.temp[paste0(var1, ".RLDC")][[1]][m1]
                    rldc2 <- dieter.temp[paste0(var2, ".RLDC")][[1]][m2]

                    # Without curtailment
                    rldc.temp <- pmax(rldc1, 0) - pmax(rldc2, 0)

                    dieter.temp <- dieter.temp %>%
                        mutate(value = rldc.temp) %>%
                        mutate(hour.sorted = 1:8760) %>%
                        mutate(technology = var2) %>%
                        select(tall, hour.sorted, technology, value) %>%
                        mutate(iteration = dieter.iter.step * i)

                    # Calculate curtailment for solar and wind
                    if (v %in% c(1, 2)) {
                        rldc.temp.curt <- pmin(rldc1, 0) - pmin(rldc2, 0)

                        dieter.temp.curt <- dieter.temp %>%
                            mutate(value = -rldc.temp.curt) %>%
                            mutate(technology = paste0(var2, "_curt"))

                        dieter.temp <- rbind(dieter.temp, dieter.temp.curt)
                    }

                    # Append output
                    out.dieter.rldc <- rbind(out.dieter.rldc, dieter.temp)
                }
            }
        }

        # Plotting ----------------------------------------------------------------

        color.mapping.rldc <- c(color.mapping, "Wind_curt" = "#66b3ff", "Solar_curt" = "#ffd940")

        swlatex(sw, "\\onecolumn")
        swlatex(sw, paste0("\\section{Residual load duration curves (RLDCs)}"))

        for (i in unique(out.dieter.rldc$iteration)) {
            plot.dieter.rldc <- out.dieter.rldc %>%
                filter(iteration == i) %>%
                mutate(technology = factor(technology, levels = c(rev(order), "Wind_curt", "Solar_curt"))) %>% # Sort technologies for plotting
                filter(hour.sorted %in% c(seq(1, 8760, 20), 8760)) %>% # Only plot every 20-th hour to decrease file size
                mutate(value = value / 1e3) # MW -> GW

            plot.dieter.demand <- out.dieter.demand %>%
                filter(iteration == i) %>%
                filter(demand.hour.sorted %in% c(seq(1, 8760, 20), 8760)) %>% # Only plot every 20-th hour to decrease file size
                mutate(demand.RLDC = demand.RLDC / 1e3) # MW -> GW

            swlatex(sw, paste0("\\subsection{RLDCs in iteration ", i, "}"))

            p <- ggplot() +
                geom_area(data = plot.dieter.rldc, aes(x = hour.sorted, y = value, fill = technology), position = "stack", alpha = 0.8) +
                scale_fill_manual(name = "Technology", values = color.mapping.rldc) +
                geom_line(data = plot.dieter.demand, aes(x = demand.hour.sorted, y = demand.RLDC, color = "Demand"), size = 1) +
                scale_color_manual(name = NULL, values = c("Demand" = "black")) +
                guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
                xlab("Hours (sorted)") +
                ylab("Generation [GW]") +
                facet_wrap(~tall, ncol = 4, scales = "free")

            swfigure(sw, print, p, sw_option = "width=20, height=12")
        }
    }

    DIETERplotPriceDurationCurves <- function() {

        # Set variables to NULL for code check compliance
        var <- tall <- value <- hour.sorted <- NULL
        model <- all_regi <- NULL

        # Data preparation --------------------------------------------------------

        cat("Plot price duration curve \n")

        # Initialise output files
        out.dieter <- NULL

        # Loop over DIETER iterations
        for (i in 1:length(dieter.files.report)) {
            # Read in demand(hour)
            dieter.report_hours <-
                file.path(outputdir, dieter.files.report[i]) %>%
                read.gdx("report_hours",
                    squeeze = F,
                    colNames = c("file", "model", "tall", "all_regi", "var", "hour", "value")
                ) %>%
                select(!c(file, model, all_regi)) %>%
                filter(var == "price") %>%
                mutate(hour = as.numeric(substring(hour, 2))) %>%
                group_by(tall) %>%
                arrange(desc(value), group_by = T) %>% # Descending order
                mutate(hour.sorted = seq(1, 8760)) %>% # Descending hour (sorted)
                mutate(iteration = dieter.iter.step * i) %>%
                ungroup()

            out.dieter <- rbind(out.dieter, dieter.report_hours)
        }

        # Plotting ----------------------------------------------------------------

        swlatex(sw, "\\onecolumn")
        swlatex(sw, paste0("\\section{Price duration curves}"))

        for (t.rep in report.periods) {
            plot.dieter <- out.dieter %>%
                filter(tall == t.rep)

            swlatex(sw, paste0("\\subsection{Price duration curve in ", t.rep, " over iterations}"))

            p <- ggplot() +
                geom_line(data = plot.dieter, aes(x = hour.sorted, y = value)) +
                xlab("Hours (sorted)") +
                ylab("Electricity price [$/MWh]") +
                scale_y_continuous(trans = "log10") +
                facet_wrap(~iteration, ncol = 4)

            swfigure(sw, print, p, sw_option = "width=20, height=10")
        }
    }

    DIETERplotInverseScreeningCurve <- function() {

        # Set variables to NULL for code check compliance
        var <- technology <- value <- tall <- hour.running <- hour.sorted <- NULL
        iteration <- screening <- inv.screening <- NULL
        IC <- OM <- FC <- CO2 <- NULL
        `CO2 cost` <- `O&M cost` <- `annualized investment cost` <- `fuel cost (divided by eta)` <- NULL
        model <- all_regi <- NULL

        # Data preparation --------------------------------------------------------

        cat("Plot inverse screening curve \n")

        dieter.report.lcoe.kW <- c("annualized investment cost", "O&M cost")
        dieter.report.lcoe.MWh <- c("fuel cost (divided by eta)", "CO2 cost")
        dieter.report.vars <- c(dieter.report.lcoe.kW, dieter.report.lcoe.MWh)

        report.tech <- c("CCGT", "OCGT", "Lignite", "Hard coal", "Biomass")

        # Initialise output files
        out.dieter <- NULL
        # Loop over DIETER iterations
        for (i in 1:length(dieter.files.report)) {
            dieter.report_tech <-
                file.path(outputdir, dieter.files.report[i]) %>%
                read.gdx("report_tech", squeeze = F, colNames = c("file", "model", "tall", "all_regi", "var", "technology", "value")) %>%
                select(!c(file, model, all_regi)) %>%
                filter(var %in% dieter.report.vars) %>%
                revalue.levels(technology = dieter.tech.mapping) %>%
                filter(technology %in% report.tech) %>%
                spread(var, value) %>%
                rename(IC = `annualized investment cost`, OM = `O&M cost`, FC = `fuel cost (divided by eta)`, CO2 = `CO2 cost`) %>%
                tidyr::replace_na(list(CO2 = 0, FC = 0))

            dieter.report_tech_hours <-
                file.path(outputdir, dieter.files.report[i]) %>%
                read.gdx("report_tech_hours",
                    squeeze = F,
                    colNames = c("file", "model", "tall", "all_regi", "var", "technology", "hour", "value")
                ) %>%
                filter(var == "generation") %>%
                select(!c(file, model, all_regi, var)) %>%
                revalue.levels(technology = dieter.tech.mapping) %>%
                mutate(hour = as.numeric(substring(hour, 2))) %>%
                filter(technology %in% report.tech) %>%
                complete(technology, tall, hour = 1:8760, fill = list(value = NA)) %>% # Fill up missing hours with 0
                group_by(tall, technology) %>%
                arrange(desc(value)) %>%
                mutate(hour.sorted = 1:8760) %>%
                mutate(hour.running = case_when(
                    is.na(value) ~ "No",
                    !is.na(value) ~ "Yes"
                )) %>%
                mutate(hour.running = factor(hour.running, levels = c("Yes", "No")))

            dieter.report <- full_join(dieter.report_tech, dieter.report_tech_hours) %>%
                mutate(screening = hour.sorted * (FC + CO2) / 1e3 + (IC + OM)) %>% # $/kW
                mutate(inv.screening = (FC + CO2) + (IC + OM) / hour.sorted * 1e3) %>% # $/MWh
                mutate(iteration = dieter.iter.step * i)

            out.dieter <- rbind(out.dieter, dieter.report)
        }

        # Plotting ----------------------------------------------------------------

        swlatex(sw, "\\onecolumn")
        swlatex(sw, paste0("\\section{Screening curves}"))

        for (i in unique(out.dieter$iteration)) {
            plot.dieter <- out.dieter %>%
                filter(iteration == i) %>%
                filter(hour.sorted %in% c(seq(1, 8760, 20), 8760)) # Only plot every 20-th hour to decrease file size

            swlatex(sw, paste0("\\subsection{Screening curves in iteration ", i, "}"))

            p <- ggplot() +
                geom_line(data = plot.dieter, aes(x = hour.sorted, y = screening, color = technology, linetype = hour.running), size = 1) +
                scale_color_manual(name = "Technology", values = color.mapping) +
                scale_linetype_manual(name = "Running", values = c(Yes = "solid", No = "dotted")) +
                xlab("Hours (sorted)") +
                ylab("Total cost [$/kW]") +
                scale_x_continuous(limits = c(0, 8760)) +
                facet_wrap(~tall)

            swfigure(sw, print, p, sw_option = "width=20, height=10")
        }

        swlatex(sw, paste0("\\section{Inverse screening curves}"))

        for (i in unique(out.dieter$iteration)) {
            plot.dieter <- out.dieter %>%
                filter(iteration == i) %>%
                filter(hour.sorted %in% c(seq(1, 8760, 20), 8760)) # Only plot every 20-th hour to decrease file size

            swlatex(sw, paste0("\\subsection{Inverse screening curves in iteration ", i, "}"))

            p <- ggplot() +
                geom_line(data = plot.dieter, aes(x = hour.sorted, y = inv.screening, color = technology, linetype = hour.running), size = 1) +
                scale_color_manual(name = "Technology", values = color.mapping) +
                scale_linetype_manual(name = "Running", values = c(Yes = "solid", No = "dotted")) +
                scale_x_continuous(limits = c(0, 8760)) +
                scale_y_continuous(limits = c(0, 200)) +
                xlab("Hours (sorted)") +
                ylab("Total cost [$/MWh]") +
                facet_wrap(~tall)

            swfigure(sw, print, p, sw_option = "width=20, height=10")
        }
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

    out <- DIETERplotCapacityFactors()
    out.dieter.capfac <- out$DIETER
    out.remind.capfac <- out$REMIND

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

    DIETERplotRLDCs()

    # Price duration curves ---------------------------------------------------

    DIETERplotPriceDurationCurves()

    # (Inverse) screening curves ----------------------------------------------

    DIETERplotInverseScreeningCurve()

    # Markups -----------------------------------------------------------------


    # Close LaTeX PDF ---------------------------------------------------------

    swclose(sw)
}