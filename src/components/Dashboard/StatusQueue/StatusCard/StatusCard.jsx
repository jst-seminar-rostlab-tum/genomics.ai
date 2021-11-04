import React from 'react';
import { Card, CardActions, CardContent, Typography, Box } from '@mui/material';

function StatusCard({ index }) {
  return (
    <Box sx={{ minWidth: '350px', backgroundColor: 'red', width: '100%'}}>
      <Card variant="outlined">
        <CardContent>
          <Typography sx={{ fontSize: 14 }} color="text.secondary" gutterBottom>
            {index}
          </Typography>
        </CardContent>
      </Card>
    </Box>
  );
}

export default StatusCard;
