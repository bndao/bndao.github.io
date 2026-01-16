---
layout: page
title: Bioinformatic analyses
permalink: /bioinformatic-analyses/
nav: true
nav_order: 2
display_categories: [Single-cell RNA-seq]
horizontal: false
---

Welcome to my GitHub portfolio dedicated to bioinformatics analyses.
This repository showcases various projects related to data analysis, statistics, and computational biology. My work focuses on reproducible research workflows and the analysis of biological data using high-performance computing resources.
Most of the computationally intensive analyses are performed using the [Mésocentre de calcul de Franche-Comté](http://meso.univ-fcomte.fr/liens.html) 
The main languages and tools I use include:
•	R for statistical analysis, data visualization, and bioinformatics packages
•	Python for data processing, scripting, and pipeline development
•	Bash for workflow automation and interaction with HPC environments
Feel free to explore the repositories, scripts, and notebooks. Feedback and discussions are welcome.

<!-- pages/projects.md -->
<div class="projects">
{% if site.enable_project_categories and page.display_categories %}
  <!-- Display categorized projects -->
  {% for category in page.display_categories %}
  <a id="{{ category }}" href=".#{{ category }}">
    <h2 class="category">{{ category }}</h2>
  </a>
  {% assign categorized_projects = site.projects | where: "category", category %}
  {% assign sorted_projects = categorized_projects | sort: "importance" %}
  <!-- Generate cards for each project -->
  {% if page.horizontal %}
  <div class="container">
    <div class="row row-cols-1 row-cols-md-2">
    {% for project in sorted_projects %}
      {% include projects_horizontal.liquid %}
    {% endfor %}
    </div>
  </div>
  {% else %}
  <div class="row row-cols-1 row-cols-md-3">
    {% for project in sorted_projects %}
      {% include projects.liquid %}
    {% endfor %}
  </div>
  {% endif %}
  {% endfor %}

{% else %}

<!-- Display projects without categories -->

{% assign sorted_projects = site.projects | sort: "importance" %}

  <!-- Generate cards for each project -->

{% if page.horizontal %}

  <div class="container">
    <div class="row row-cols-1 row-cols-md-2">
    {% for project in sorted_projects %}
      {% include projects_horizontal.liquid %}
    {% endfor %}
    </div>
  </div>
  {% else %}
  <div class="row row-cols-1 row-cols-md-3">
    {% for project in sorted_projects %}
      {% include projects.liquid %}
    {% endfor %}
  </div>
  {% endif %}
{% endif %}
</div>
