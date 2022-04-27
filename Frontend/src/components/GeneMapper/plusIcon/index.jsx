import Icon from '@mui/material/Icon';
import React from 'react';
import './plusIcon.css';

export default function PlusIcon() {
  return (
    // change the size of it
    <Icon className="plusicon" color="primary" sx={{ fontSize: 53 }}>add_circle</Icon>

  );
}
