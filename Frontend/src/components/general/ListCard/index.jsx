import React from 'react';
import Card from '@mui/material/Card';
import styles from './listCard.module.css';

function InstitutionCard({
  imageComponent, imageURL, title, description, nextToTitle, trailing,
}) {
  if (imageURL) {
    if (imageComponent) {
      throw new Error('You can only have either imageComponent or imageURL.');
    }
    // TODO: change imageComponent
  }
  return (
    <Card sx={{ background: 'card.background' }}>
      <div className={styles.cardContent}>
        <div className={styles.start}>
          {imageURL && (
            <img className={styles.cardImage} src={imageURL} alt="" />
          )}
          <div className={styles.text}>
            <div className={styles.titleRow}>
              <h3 className={styles.title}>{title}</h3>
              {nextToTitle || (<></>)}
            </div>
            <div className={styles.description}>
              {description && (
                <p>{description}</p>
              )}
            </div>
          </div>
        </div>
        {trailing && (
          <div className={styles.trailing}>
            {trailing}
          </div>
        )}
      </div>
    </Card>
  );
}

export default InstitutionCard;
