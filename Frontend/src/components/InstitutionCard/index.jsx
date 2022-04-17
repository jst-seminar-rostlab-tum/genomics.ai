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
        <div>
          <img src={profilePictureURL} alt="" />
          <div>
            <h3>{name}</h3>
            <p>{country}</p>
          </div>
        </div>
      </CardContent>
    </Card>
  );
}

export default InstitutionCard;
