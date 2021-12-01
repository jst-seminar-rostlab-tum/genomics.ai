import React from 'react';
import { Typography, Box } from '@mui/material';
import NavBar from '../../../NavBar/NavBar';
import styles from './home.module.css';
import logo1 from '../../../../assets/landing-illustrations/science.png';
import logo2 from '../../../../assets/landing-illustrations/files.png';
import logo3 from '../../../../assets/landing-illustrations/dashboard.png';
import logo4 from '../../../../assets/landing-illustrations/visual.png';
import Footer from '../../Footer/Footer';
import dnaImage from '../../../../assets/dna.png';

function Home(props) {
  const { setUser } = props;
  return (
    <div className={styles.container}>
      <NavBar setUser={setUser} />
      <img className={styles.backgroundImage} src={dnaImage} alt="background" />
      {/* If genomics visualized stays in the middle, then it needs to be properly aligned.
      Right now it is using padding as the way to align the items */}
      <Box sx={{ paddingBottom: '50px', paddingTop: '250px' }}>
        <Typography sx={{ fontSize: '24px', fontWeight: 'bold' }}>Genomics visualized</Typography>
      </Box>
      <div className={styles.infoContainer}>
        <img className={styles.illustration1} src={logo1} alt="science-guy" />
        <div className={styles.explanationRight}>
          <Typography sx={{ fontSize: '24px' }}>Visualize all of your single-cell sequencing data</Typography>
          <Typography sx={{ fontSize: '24px' }}>in a fast and easy way</Typography>
        </div>
      </div>
      <div className={styles.infoContainer}>
        <div className={styles.explanationLeft}>
          <Typography className={styles.illustrationTitle} sx={{ fontSize: '24px', fontWeight: 'bold' }}>Upload</Typography>
          <Typography sx={{ fontSize: '24px' }}>Upload your files containing the single-cell</Typography>
          <Typography sx={{ fontSize: '24px' }}>sequencing data</Typography>
        </div>
        <img className={styles.illustration2} src={logo2} alt="explanation" />
      </div>
      <div className={styles.infoContainer}>
        <img className={styles.illustration3} src={logo3} alt="explanation" />
        <div className={styles.explanationRight}>
          <Typography className={styles.illustrationTitle} sx={{ fontSize: '24px', fontWeight: 'bold' }}>Create a project</Typography>
          <Typography sx={{ fontSize: '24px' }}>Create and access all of your projects</Typography>
          <Typography sx={{ fontSize: '24px' }}>in the dashboard</Typography>
        </div>
      </div>
      <div className={styles.infoContainer}>
        <div className={styles.explanationLeft}>
          <Typography className={styles.illustrationTitle} sx={{ fontSize: '24px', fontWeight: 'bold' }}>See results</Typography>
          <Typography sx={{ fontSize: '24px' }}>See your results after the algorithm</Typography>
          <Typography sx={{ fontSize: '24px' }}>has processed the data</Typography>
        </div>
        <img className={styles.illustration4} src={logo4} alt="explanation" />
      </div>
      <Footer />
    </div>
    // TODO: add footer for the website across all pages that are not the dashboard

  );
}

export default Home;
