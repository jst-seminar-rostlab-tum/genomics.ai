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

function loadingOrErrorOrResult(object) {
  if (!object) {
    return <Typography sx={{ paddingTop: '10rem' }} fontSize="5rem">Fetching result...</Typography>;
  } if (object.error) {
    return <Typography sx={{ paddingTop: '10rem' }} fontSize="5rem">{object.error}</Typography>;
  }
  return <Visualization />;
}

function VisualizationPage() {
  const [job, setJob] = useState(null);
  const [jobId, setJobId] = useState(null);
  const query = useQuery();

  useEffect(() => {
    const id = query.get('id');
    setJobId(id);
    queryJob(id).then((response) => setJob(response))
      .catch((err) => {
        console.log(err);
        console.log('WTF');
        setJob({ error: 'Couldn\'t fetch job result.' });
      });
  }, [query]);
  return (
    <ThemeProvider theme={theme}>
      <div>
        <NavBar />
        <div className={styles.headerContainer}>
          <Typography sx={{ paddingTop: '8rem' }} fontSize="3rem">
            Job result
            <br />
            {jobId}
          </Typography>
          <Box
            sx={{
              display: 'flex',
              justifyContent: 'center',
              textAlign: 'center',
            }}
            className={styles.visualizationContainer}
          >
            {
              loadingOrErrorOrResult(job)
            }
          </Box>
          <Footer />
        </div>
      </div>
    </ThemeProvider>
  );
}

export default VisualizationPage;
