# Base R Shiny image
FROM mziemann/gmea

# Update codes
RUN cd /gmea && git pull

# Make a directory in the container
RUN mkdir /home/app

# Copy the Shiny app code
RUN cp -r /gmea/app /home/

# Expose the application port
EXPOSE 3838

# Set the container working directory
ENV DIRPATH /home/app
WORKDIR $DIRPATH

# Run the R Shiny app
CMD Rscript /home/app/app.R
