import React from 'react';
import {
  Box, createTheme, ThemeProvider,
} from '@mui/material';
import { useLocation } from 'react-router-dom';
import Visualization from '../../components/Visualization/src/Visualization';
import styles from './visualizationpage.module.css';
import NavBar from '../../components/NavBar';
import Footer from '../../components/Footer';

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
