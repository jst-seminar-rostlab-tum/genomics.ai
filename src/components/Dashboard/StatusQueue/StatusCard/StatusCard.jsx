import React, { useState } from 'react';
import { createTheme, ThemeProvider } from '@mui/material/styles';
import {
  Box,
  Typography,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Button,
  Stack,
} from '@mui/material';
import {
  red,
  green,
  yellow,
  grey,
  blue,
} from '@mui/material/colors';
import CircleIcon from '@mui/icons-material/Circle';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import styles from './statuscard.module.css';

/*
The status card includes the status of the job.
Status types that the status card can include are:

1. pending - state during the upload of the file and during verification

2. processing - state during the processing of the file for the ML model

3. error - can occure in all steps including verification

4. completed - state after all stages are completed

5. unknown - the unknown state happens when the status can't be updated
*/
const theme = createTheme({
  palette: {
    primary: {
      main: '#1888ff',
    },
    secondary: {
      main: '#6B7379',
    },
  },
});

function StatusCard({ id }) {
  const [response, setResponse] = useState({
    status: 'unknown',
    log: "The file's status is currently unknown.",
  });

  // object of colors for the statuses
  const statusColor = {
    pending: yellow[600],
    processing: blue[300],
    error: red[300],
    completed: green[300],
    unknown: grey[500],
  };

  // Change of state for testing
  // The status is not passed down as a prop, but received from backend,
  const [changeResCount, setChangeResCount] = useState(0);
  if (changeResCount < 1) {
    setResponse({ ...response, status: 'error' });
    setChangeResCount(changeResCount + 1);
  }

  return (
    <ThemeProvider theme={theme}>
      <Box className={styles.cardContainer}>
        <Accordion>
          <AccordionSummary
            expandIcon={<ExpandMoreIcon />}
            aria-controls="panel1a-content"
            id="panel1a-header"
          >
            <div className={styles.summaryContainer}>
              <Stack
                spacing={2}
                direction="row"
              >
                <Typography
                  className={styles.headerContainer}
                  variant="h4"
                  sx={{
                    fontWeight: 'bold',
                    fontSize: '16px',
                    margin: '0',
                    padding: '0',
                    textAlign: 'center',
                  }}
                >
                  {`Job ${id}`}
                </Typography>
                <CircleIcon sx={{ color: statusColor[response.status], height: '20px' }} />
                <Typography
                  sx={{
                    fontSize: '16px',
                    fontWeight: 'light',
                    margin: '0',
                    padding: '0',
                    textAlign: 'center',
                    display: 'inline',
                    color: statusColor[response.status],
                  }}
                >
                  {` ${response.status.toUpperCase()}`}
                </Typography>
              </Stack>

            </div>
          </AccordionSummary>
          <AccordionDetails className={styles.fileStatusLog}>
            <div className={styles.statusContainer}>
              <Typography variant="h4" sx={{ fontWeight: 'light', fontSize: '16px' }}>
                Status:
                <Typography
                  sx={{
                    fontSize: '16px',
                    fontWeight: 'light',
                    margin: '0',
                    padding: '0',
                    textAlign: 'center',
                    display: 'inline',
                    color: statusColor[response.status],
                  }}
                >
                  {` ${response.status.toUpperCase()}`}
                </Typography>
              </Typography>
            </div>
            <Typography>{response.log}</Typography>

            <Box sx={{ textAlign: 'left', paddingTop: '30px' }}>
              { response.status === 'completed' ? (
                <Stack
                  spacing={2}
                  direction="row"
                >
                  <Button
                    href="#"
                    variant="outlined"
                    color="secondary"
                    style={{ borderRadius: '10px' }}
                  >
                    Delete
                  </Button>

                  <Button
                    href="#"
                    variant="contained"
                    color="primary"
                    style={{ borderRadius: '10px' }}
                  >
                    View Result
                  </Button>

                </Stack>
              )
                : (
                  <Button
                    href="#"
                    variant="outlined"
                    color="secondary"
                    style={{ borderRadius: '10px' }}
                  >
                    Cancel
                  </Button>
                )}
            </Box>

          </AccordionDetails>
        </Accordion>
      </Box>
    </ThemeProvider>

  );
}

export default StatusCard;
