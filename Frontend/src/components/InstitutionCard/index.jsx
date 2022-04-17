import React from 'react';
import Card from '@mui/material/Card';
import CardContent from '@mui/material/CardContent';

function InstitutionCard({ institution }) {
  const {
    id, name, country, profilePictureURL, backgroundPictureURL, adminIds,
  } = institution;
  return (
    <Card raised>
      <CardContent>
        <h2>{name}</h2>
      </CardContent>
    </Card>
  );
}

export default InstitutionCard;
