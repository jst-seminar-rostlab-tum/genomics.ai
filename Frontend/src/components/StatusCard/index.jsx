import React from 'react';
import { Link } from 'react-router-dom';
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
import CircleIcon from '@mui/icons-material/Circle';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import { jobStatusColors, jobStatusTitles } from 'shared/utils/common/constants';
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

function StatusCard({
  id, status, log, location,
}) {
  // TODO: talk to Dom to check with prepending the path to the results page is necessary

  const cancelJob = () => {
    // TODO: cancel the job here
  };

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
                  {`Job ${id.substr(id.length - 4)}`}
                </Typography>
                <CircleIcon sx={{ color: jobStatusColors[status], height: '20px' }} />
                <Typography
                  sx={{
                    fontSize: '16px',
                    fontWeight: 'bold',
                    margin: '0',
                    padding: '0',
                    textAlign: 'center',
                    display: 'inline',
                    color: jobStatusColors[status],
                  }}
                >
                  {` ${jobStatusTitles[status]}`}
                </Typography>
              </Stack>

            </div>
          </AccordionSummary>
          <AccordionDetails className={styles.fileStatusLog}>
            <Box sx={{ textAlign: 'center' }}>
              <Typography variant="h4" sx={{ fontWeight: '500', fontSize: '16px' }}>
                Status:
                <Typography
                  sx={{
                    fontSize: '16px',
                    fontWeight: 'bold',
                    margin: '0',
                    padding: '0',
                    textAlign: 'center',
                    display: 'inline',
                    color: jobStatusColors[status],
                  }}
                >
                  {` ${jobStatusTitles[status]}`}
                </Typography>
              </Typography>
              <Typography sx={{ paddingBottom: '10px' }}>
                id:
                {' '}
                {id}
              </Typography>
              <Typography>{log}</Typography>
            </Box>
            {/* Results button
            //TODO: the href="" in the button should point to the right visualization result
            setting currently to "result" + id number, i.e. : "result11s2ef2"
            */}
            <Box sx={{ textAlign: 'center', paddingTop: '30px' }}>
              {status === 'DONE' ? (
                <Link
                  to={{ pathname: '/result', search: `tsv=${location}` }}
                  target="_blank"
                  rel="noopener noreferrer"
                  style={{ textDecoration: 'none' }}
                >
                  <Button variant="outlined" color="info">
                    {/* The id is given as a prop */}
                    <Typography sx={{ color: '#4F83CC', fontWeight: '500' }}>
                      See results
                    </Typography>
                  </Button>
                </Link>
              ) : null}
              {/* cancel button */}
              {status !== 'ABORTED' && status !== 'DONE' && (
                <Button variant="outlined" color="error" onClick={cancelJob}>
                  cancel
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
