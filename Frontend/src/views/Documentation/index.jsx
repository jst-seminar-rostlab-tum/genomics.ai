import React from 'react';
import { Divider, Stack, Typography } from '@mui/material';
import styles from './documentation.module.css';
import visualisation from '../../assets/landing-illustrations/visualisation.png';
import csvformat from '../../assets/landing-illustrations/csvformat.png';
import zoom from '../../assets/landing-illustrations/zoomVisualisation.png';
import clusterLabel from '../../assets/landing-illustrations/clusterLabel.png';
import allClusterLabels from '../../assets/landing-illustrations/allClusterLabels.png';
import clusterDetails from '../../assets/landing-illustrations/clusterDetails.png';

function Documentation() {
  return (
    <div className={styles.container}>
      <Stack
        spacing="30px"
      >
        {/* Documentation Headline ---------------------------------------------------------*/}
        <Typography sx={{ fontWeight: 'bold', fontSize: '30px' }}>Documentation</Typography>
        <div className={styles.textContainer}>
          <Typography
            sx={{ fontSize: '25px', paddingInline: '350px', maxWidth: '1700px' }}
            align="center"
          >
            This Documentationn provides you with a basic understanding of our
            implemented solution.
          </Typography>
        </div>

        {/* Machine Learning Model-----------------------------------------------------------*/}
        <Divider
          variant="middle"
          textAlign="center"
          sx={{ paddingTop: '100px', paddingBottom: '40px' }}
        >
          <Typography sx={{ width: '400px', fontWeight: 'bold', fontSize: '30px' }}>
            1. Machine Learning Model
          </Typography>
        </Divider>
        <div className={styles.textContainer}>
          <Typography
            sx={{ fontSize: '20px', paddingInline: '350px', maxWidth: '1700px' }}
            align="center"
          >
            We are working hard on providing you with a documentation for our
            machine learning model in the near future.
            If you experience problems or have any questions
            feel free to get in touch with us.
          </Typography>
        </div>

        {/* Visualisation -------------------------------------------------------------------*/}
        <Divider
          variant="middle"
          textAlign="center"
          sx={{ paddingTop: '100px', paddingBottom: '40px' }}
        >
          <Typography sx={{ width: '250px', fontWeight: 'bold', fontSize: '30px' }}>2. Visualisation</Typography>
        </Divider>
        <div className={styles.textContainer}>
          <Typography
            sx={{ fontSize: '20px', paddingInline: '350px', maxWidth: '1700px' }}
            align="center"
          >
            Based on your input data set we will provide you with a visualisation containing
            the calculated classifications for your cells. (An example can be seen below)
          </Typography>
        </div>
        <img className={styles.mainIllustration} src={visualisation} alt="visualisation" />

        {/* Input Data -----------------------------------------------------------------------*/}
        <Divider
          variant="middle"
          textAlign="center"
          sx={{ paddingTop: '100px', paddingBottom: '40px' }}
        >
          <Typography sx={{ width: '150px', fontWeight: 'bold', fontSize: '30px' }}>
            Input Data
          </Typography>
        </Divider>
        <div className={styles.textContainer}>
          <Typography
            sx={{ fontSize: '20px', paddingInline: '350px', maxWidth: '1700px' }}
            align="center"
          >
            After processing the raw data provided by the user,
            the machine learning model returns a .csv file.
            This .csv file is needed as input for the visualisation in the following format

          </Typography>
        </div>

        <div className={styles.textContainer}>
          <Typography
            sx={{
              fontSize: '20px', paddingInline: '350px', maxWidth: '1700px', fontWeight: 'bold',
            }}
            align="center"
          >
            (x,y denote the UMAP coordinates of a sample)
          </Typography>
        </div>
        <img className={styles.csvImage} src={csvformat} alt="visualisation" />

        {/* Features ---------------------------------------------------------------------------*/}
        <Divider
          variant="middle"
          textAlign="center"
          sx={{ paddingTop: '100px', paddingBottom: '40px' }}
        >
          <Typography sx={{ fontWeight: 'bold', fontSize: '30px' }}>
            Features
          </Typography>
        </Divider>

        <Stack
          direction="column"
          spacing={25}
        >
          {/* Zoom in/out -----------------------------------------------------------------------*/}

          <div className={styles.infoContainer}>
            <img className={styles.verticalIllustration} src={zoom} alt="science-guy" />
            <div className={styles.explanationLeft}>
              <Typography
                className={styles.explanationTitle}
                sx={{ fontSize: '25px', fontWeight: 'bold' }}
              >
                Zoom in/out
              </Typography>

              <Typography sx={{ fontSize: '20px' }}>
                The buttons in the top right corner of the visualisation
                allow you to zoom in on the generated sample points.
              </Typography>
              <br />
              <Typography sx={{ fontSize: '20px' }}>
                At the highest zoom level, the sample IDs are also displayed when hovering.
                Panning is only possible from the first zoom level.
              </Typography>
              <br />
              <br />
              <Typography
                className={styles.explanationTitle}
                sx={{ fontSize: '25px', fontWeight: 'bold' }}
              >
                Resetting the Zoom Level
              </Typography>

              <Typography sx={{ fontSize: '20px' }}>
                By clicking the &apos;Reset&apos; button in the upper right corner
                you can reset the zoom level to the original value.
              </Typography>

            </div>
          </div>

          {/* Cluster Labels ------------------------------------------------------------------*/}
          <Stack direction="column">
            <div className={styles.infoContainer}>
              <img
                className={styles.horizontalIllustration}
                src={clusterLabel}
                alt="cluster label"
              />
              <div className={styles.explanationLeft}>
                <Typography
                  className={styles.explanationTitle}
                  sx={{ fontSize: '25px', fontWeight: 'bold' }}
                >
                  Single Cluster Labels
                </Typography>
                <Typography sx={{ fontSize: '20px' }}>
                  {' '}
                  The visualization temporarily displays a cluster’s name when
                  hovering over it.
                </Typography>
              </div>
            </div>

            <div className={styles.infoContainer}>
              <img
                className={styles.horizontalIllustration}
                src={allClusterLabels}
                alt="all cluster labels"
              />
              <div className={styles.explanationLeft}>
                <Typography
                  className={styles.explanationTitle}
                  sx={{ fontSize: '25px', fontWeight: 'bold' }}
                >
                  All Cluster Labels
                </Typography>
                <Typography sx={{ fontSize: '20px' }}>
                  {' '}
                  To permanently display all cluster labels you need to press the &apos;Labels&apos;
                  button in the right upper corner.
                </Typography>
              </div>
            </div>
          </Stack>

          {/* Cluster Details ------------------------------------------------------------------*/}
          <div className={styles.infoContainer}>
            <img
              className={styles.horizontalIllustration}
              src={clusterDetails}
              alt="cluster Details"
            />
            <div className={styles.explanationLeft}>
              <Typography
                className={styles.explanationTitle}
                sx={{ fontSize: '25px', fontWeight: 'bold' }}
              >
                Cluster Details
              </Typography>

              <Typography sx={{ fontSize: '20px' }}>
                When clicking on a cluster you’ll be presented with more
                detailed information about it.
              </Typography>
              <br />
              <br />
              <Typography
                className={styles.explanationTitle}
                sx={{ fontSize: '25px', fontWeight: 'bold' }}
              >
                Portrayed Cluster Attributes:
              </Typography>
              <Typography sx={{ fontSize: '20px' }}>
                - percentage share of all samples
              </Typography>
              <Typography sx={{ fontSize: '20px' }}>
                - number of samples
              </Typography>
              <Typography sx={{ fontSize: '20px' }}>
                - average L2 distance of the samples to the mean.
              </Typography>

            </div>
          </div>
        </Stack>

      </Stack>
    </div>
  );
}

export default Documentation;
