import React, { useState } from 'react';
import {
  Card,
  Box,
  Typography,
  Collapse,
  IconButton,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Button
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
import ExpandLessIcon from '@mui/icons-material/ExpandLess';
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

function StatusCard({ id }) {
  const [response, setResponse] = useState({
    status: 'unknown',
    log: "The file's status is currently unknown.",
  })

  const [moreInfo, setMoreInfo] = useState(false);

  // object of colors for the statuses
  const statusColor = {
    pending: yellow[600],
    processing: blue[300],
    error: red[300],
    completed: green[300],
    unknown: grey[500],
  };

  // The status is not passed down as a prop, but received from backend
  const [changeResCount, setChangeResCount] = useState(0);
  if (changeResCount < 1) {
    setResponse({ ...response, status: 'completed' });
    setChangeResCount(changeResCount + 1);
  } // check dom's code to understand how it is done

  return (
    <Box className={styles.cardContainer}>

      <Accordion>
        <AccordionSummary
          expandIcon={<ExpandMoreIcon />}
          aria-controls="panel1a-content"
          id="panel1a-header"
        >
          <div className={styles.summaryContainer}>
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
            <CircleIcon sx={{ color: statusColor[response.status], paddingLeft: '5px' }} />
          </div>
        </AccordionSummary>
        <AccordionDetails className={styles.fileStatusLog}>
          <div className={styles.statusContainer}>
            <Typography variant="h4" sx={{ fontWeight: 'light', fontSize: '16px' }}>
              Status:
              <Typography
                variant="h4"
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
          <Box sx={{ textAlign: 'center', paddingTop: '20px' }}>
            { response.status === 'completed' ? <Button href="#">View results</Button> : null }
          </Box>
        </AccordionDetails>
      </Accordion>
    </Box>
  );
}

export default StatusCard;

//using accordion instead of current implementation

/*
<div className="infoContainer">
          <Typography variant="h4" sx={{ fontWeight: 'light', fontSize: '16px' }}>
            Status:
            <Typography
              variant="h4"
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
          <Typography>hello</Typography>
          <IconButton onClick={() => setMoreInfo(!moreInfo)}>
            {moreInfo ? <ExpandLessIcon /> : <ExpandMoreIcon />}
          </IconButton>
        </div>
*/


/*
<Accordion>
          <AccordionSummary
            expandIcon={<ExpandMoreIcon />}
            aria-controls="panel1a-content"
            id="panel1a-header"
          >
            <Typography>Accordion 1</Typography>
          </AccordionSummary>
          <AccordionDetails>
            <Typography>
              Lorem ipsum dolor sit amet, consectetur adipiscing elit. Suspendisse
              malesuada lacus ex, sit amet blandit leo lobortis eget.
            </Typography>
          </AccordionDetails>
        </Accordion>
*/

/**
 *
 * <Card variant="outlined" sx={{ padding: '10px' }}>
        <div className={styles.headerContainer}>
          <Typography
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
          <CircleIcon sx={{ color: statusColor[response.status], paddingLeft: '5px' }} />
        </div>
        <LogAccordion expand={moreInfo} />
        <Accordion>
          <AccordionSummary
            expandIcon={<ExpandMoreIcon />}
            aria-controls="panel1a-content"
            id="panel1a-header"
          >
            <Typography>Accordion 1</Typography>
          </AccordionSummary>
          <AccordionDetails>
            <Typography>
              Lorem ipsum dolor sit amet, consectetur adipiscing elit. Suspendisse
              malesuada lacus ex, sit amet blandit leo lobortis eget.
            </Typography>
          </AccordionDetails>
        </Accordion>
      </Card>
 * / */