import { Typography, Box } from '@mui/material';
import React from 'react';

// Component to display the feedback text
const ResultStatus = ({ count, searchedEntity, searchedKeyword }) => {
  const countText = `${count} result${count === 1 ? '' : 's'}`;
  const entityText = searchedEntity.charAt(0).toUpperCase() + searchedEntity.slice(1);
  const keywordText = searchedKeyword && searchedKeyword !== '' && (
    <>
      matching
      {' '}
      <b>{searchedKeyword}</b>
    </>
  );
  return (
    <Box sx={{ m: 2 }}>
      <Typography>
        {`${countText} for ${entityText} `}
        {keywordText}
      </Typography>
    </Box>
  );
};

export default ResultStatus;
