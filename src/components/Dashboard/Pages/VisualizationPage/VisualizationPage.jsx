import React, { useEffect, useState } from 'react';
import {
  Box, createTheme, ThemeProvider, Typography,
} from '@mui/material';
import { useLocation } from 'react-router-dom';
import Visualization from '../../Visualization/src/Visualization';
import styles from './visualizationpage.module.css';
import queryJob from './VisualizationPageLogic';
import NavBar from '../../../NavBar/NavBar';
import Footer from '../../../LandingPage/Footer/Footer';

const theme = createTheme({
  palette: {
    primary: {
      main: '#0075FF',
    },
    secondary: {
      main: '#FFFFFF',
    },
  },
});

function useQuery() {
  const { search } = useLocation();

  return React.useMemo(() => new URLSearchParams(search), [search]);
}

function VisualizationPage() {
  const [job, setJob] = useState(null);
  const [jobId, setJobId] = useState(null);
  const query = useQuery();
  return (
    <ThemeProvider theme={theme}>
      <div>
        <NavBar />
        <div className={styles.headerContainer}>
          <Box
            sx={{
              display: 'flex',
              justifyContent: 'center',
              textAlign: 'center',
            }}
            className={styles.visualizationContainer}
          >
            <Visualization url={query.get('tsv')} />
          </Box>
          <Footer />
        </div>
      </div>
    </ThemeProvider>
  );
}

export default VisualizationPage;
