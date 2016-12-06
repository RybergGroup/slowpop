#' A slowpop function
#'
#' This function creates a population
#' @param n number of individuals in the population
#' @param diploid is the organism diploid or not (TRUE/FALSE)
#' @param haplotypes the number of haplotypes in the population
#' @keywords population
#' @export
#' @examples
#' my_population <- create_population(100)

create_population <- function(n,diploid=TRUE,haplotypes=2) {
    population<-list()
    for (i in 1:n) {
	if (diploid) population[[i]]<-(sample(1:haplotypes,2,replace=T)) else
	    population[[i]]<-(sample(1:haplotypes,1))
    }
    return(population)
}

#' A slowpop function
#'
#' This function test if the population exist
#' @param population is the population to test
#' @keywords population test
#' @export
#' @examples
#' my_population <- create_population(100)
#' test_population(my_population)

test_population<-function(population) {
    if (!is.list(population) || length(population)<1) return(FALSE) else
	return(TRUE)
}

#' A slowpop function
#'
#' Calculates the relative probability of each individual to make it into the next generation
#' @param population is the population
#' @param relative_fitness a named vector with the relative fitness of each haplotype
#' @keywords population test
#' @export
#' @examples
#' my_population <- create_population(100)
#' my_relative_fitness <- get_named_vector(2)
#' my_relative_fitness <- c(1,1.1,1.1)
#' my_probabilities <- assign_relative_fitness(my_population,my_relative_fitness)

assign_relative_fitness<-function(population,relative_fitness){
    if (!test_population(population)) return(NA)
    probabilities<-c()
    for(i in 1:length(population)) {
	rel_f<-relative_fitness[paste(sort(population[[i]]),collapse=',')]
	#print(rel_f)
	if(length(rel_f) <2 && is.na(rel_f)) return(NA)
	probabilities[i]<-rel_f
    }
    return(probabilities)
}

#' A slowpop function
#'
#' Select individuals based on relative fitness
#' @param population is the population
#' @param relative_fitness is a vector with the relative fitness of each individual
#' @keywords population test
#' @export
#' @examples
#' my_population <- create_population(100)
#' my_relative_fitness <- get_named_vector(2)
#' my_relative_fitness <- c(1,1.1,1.1)
#' my_probabilities <- assign_relative_fitness(my_population,my_relative_fitness)


get_individuals<-function(relative_fitness,n=1) {
    if(length(relative_fitness) <2 && is.na(relative_fitness)) return(NA)
    choosen<-c()
    probabilities<-relative_fitness/sum(relative_fitness)
    #if (sum(probabilities) < 1.00) return(NA)
    for (i in 1:n) {
	number<-runif(1,0,1)
	for (j in 1:length(probabilities)) {
	    number<-number-probabilities[j]
	    choosen[i]<-j
	    if (is.na(number)) {
		cat('Number error!\n')
		if(is.na(probabilities[j])) cat('Error in calc of probabilities\n')
		return(choosen)
	    }
	    if (number <= 0) {
		probabilities[j]<-0
		probabilities<-probabilities/sum(probabilities)
		break
	    }
	}
    }
    if(length(choosen)>1 && choosen[1]==choosen[2]) cat('Error! Non outcrossing\n')
    return(choosen)
}

#' A slowpop function
#'
#' Generates a new generation from a given population
#' @param population is the population
#' @param relative_fitness a named vector with the relative fitness of each haplotype
#' @keywords population
#' @export
#' @examples
#' my_population <- create_population(100)
#' my_relative_fitness <- get_named_vector(2)
#' my_relative_fitness <- c(1,1.1,1.1)
#' my_new_generation <- new_generation(my_population,my_relative_fitness)

new_generation<-function(population,relative_fitness=NA) {
    if (!test_population(population)) {
	cat('Inappropriate population.\n')
	return(NA)
    }
    if (length(relative_fitness) <2 && is.na(relative_fitness)) {
	relative_fitness<-genotype_numbers(population)
	for(i in 1:length(relative_fitness)) relative_fitness[i]<-1
    }
    #print(relative_fitness)
    rel_fit_indv<-assign_relative_fitness(population,relative_fitness)
    if (length(rel_fit_indv)<2 && is.na(rel_fit_indv)) {
	cat('Failed to assign relative fitness to individuals.\n')
	return(NA)
    }
    #print(rel_fit_indv)
    diploid<-c()
    if (length(population[[1]])>1) diploid<-T else
	diploid<-F
    new_population<-list()
    for(i in 1:length(population)) {
	parents<-get_individuals(rel_fit_indv,2)
	if (length(parents)<2 && is.na(parents)) {
	    cat('Failed to identify parents.\n')
	    return(NA)
	}
    	new_population[[i]]<-sample(population[[parents[1]]],1)
	if(diploid) {
	    new_population[[i]][2]<-sample(population[[parents[2]]],1)
	}
    }
    return(new_population)
}

#' A slowpop function
#'
#' Get the number of homozygotes in the population
#' @param population is the population
#' @keywords homozygotes
#' @export
#' @examples
#' my_population <- create_population(100)
#' n_homozygots(my_population)

n_homozygots<-function(population) {
    if (!test_population(population)) return(NA)
    if (!(length(population[[1]])>1)) return(NA)
    n<-0
    for(i in 1:length(population)) {
	if (length(unique(population[[i]]))<2) n<-n+1
    }
    return(n)
}

#' A slowpop function
#'
#' Get the number of haplotypes in the population
#' @param population is the population
#' @keywords homozygotes
#' @export
#' @examples
#' my_population <- create_population(100)
#' n_haplotypes(my_population)

n_haplotypes<-function(population) {
    if (!test_population(population)) return(NA)
    uniques <- unique(population[[1]])
    if(length(population)>1)
	for(i in 2:length(population)) {
	    uniques<-c(uniques,population[[i]][!(population[[i]] %in% uniques)])
	}
    return(length(uniques))
}

#' A slowpop function
#'
#' Get the highest number assigned to a haplotype in the population
#' @param population is the population
#' @keywords homozygotes
#' @export
#' @examples
#' my_population <- create_population(100)
#' max_haplotype_no(my_population)

max_haplotype_no<-function(population) {
    if(!test_population(population)) return(NA)
    max_value<-0
    for(i in 1:length(population)) {
	max_value<-max(c(population[[i]],max_value))
    }
    return(max_value)
}

#' A slowpop function
#'
#' Introduce a new haplotype in the population by changing one haplotype in one individual. The haplotype get the maximum number of any existing haplotype plus one.
#' @param population is the population
#' @keywords homozygotes
#' @export
#' @examples
#' my_population <- create_population(100)
#' introduce_mutation(my_population)

introduce_mutation<-function(population) {
    if(!test_population(population)) return(NA)
    mutant<-sample(1:length(population),1)
    population[[mutant]][sample(1:length(population[[mutant]]),1)]<-max_haplotype_no(population)+1
    return(population)
}

#' A slowpop function
#'
#' Get the number of individuals for each genotype
#' @param population is the population
#' @keywords haplotypes genotypes
#' @export
#' @examples
#' my_population <- create_population(100)
#' genotype_numbers(my_population)

genotype_numbers<-function(population) {
    if(!test_population(population)) return(NA)
    genotypes<-vector()
    for (i in 1:length(population)) {
	name<-paste(sort(population[[i]]),collapse=',')
	if (name %in% names(genotypes)) {
	    genotypes[name]<-genotypes[name]+1
       	} else {
	    temp <- c(1)
	    names(temp)<- name
	    genotypes<-c(genotypes,temp)
      	}
    }
    return(genotypes)
}

#' A slowpop function
#'
#' Get the how many copies there are of each haplotype
#' @param population is the population
#' @keywords haplotypes alleles genotypes
#' @export
#' @examples
#' my_population <- create_population(100)
#' allele_numbers(my_population)

allele_numbers<-function(population,n_alleles=NA) {
    if(!test_population(population)) return(NA)
    alleles<-c()
    for(i in 1:length(population)) {
	for (j in 1:length(population[[i]])) {
	    if(is.null(alleles[population[[i]][j]]) || is.na(alleles[population[[i]][j]])) alleles[population[[i]][j]]<-1 else
		alleles[population[[i]][j]]<-alleles[population[[i]][j]]+1
	}
    }
    length<-c()
    if (is.na(n_alleles)) length<-length(alleles) else
	length<-n_alleles
    for(i in 1:length) if(is.null(alleles[i]) || is.na(alleles[i])) alleles[i]<-0
    return(alleles)
}

#' A slowpop function
#'
#' Get a vector for each combination of n haplotypes
#' @param n is the number of haplotypes
#' @keywords haplotypes alleles genotypes
#' @export
#' @examples
#' my_population <- create_population(100)
#' get_named_vector(2)

get_named_vector<-function(n=NA){
    if (is.na(n)) return(NA)
    names<-vector()
    return_vector<-vector()
    for (i in 1:n)
	for(j in i:n) {
	    names<-c(names,paste(c(i,j),collapse=','))
	    return_vector<-c(return_vector,0)
	}
    names(return_vector)<-names
    return(return_vector)
}
